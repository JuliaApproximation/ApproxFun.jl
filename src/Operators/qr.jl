



mutable struct QROperator{CO,MT,T} <: Operator{T}
    R_cache::CO
    H::MT # Contains the Householder reflections
    ncols::Int   # number of cols already upper triangularized
end

QROperator(R::CachedOperator,H::AbstractArray,ncs::Int) =
    QROperator{typeof(R),typeof(H),eltype(H)}(R,H,ncs)


convert(::Type{Operator{T}},QR::QROperator) where {T} =
    QROperator(convert(Operator{T},QR.R_cache), convert(AbstractArray{T}, QR.H),QR.ncols)

qr(QR::QROperator) = QR
factorize(QR::QROperator) = QR

for OP in (:domainspace,:rangespace)
    @eval $OP(QR::QROperator) = $OP(QR.R_cache)
end

getindex(QR::QROperator, k::Integer, j::Integer) = QR.R_cache.op[k,j]



struct QROperatorR{QRT,T} <: Operator{T}
    QR::QRT
end

QROperatorR(QR) = QROperatorR{typeof(QR),eltype(QR)}(QR)
domainspace(R::QROperatorR) = domainspace(R.QR)
rangespace(R::QROperatorR) = ℓ⁰

function getindex(R::QROperatorR,k::Integer,j::Integer)
    if j < k
        zero(eltype(R))
    else
        resizedata!(R.QR,:,j)
        R.QR.R_cache[k,j]
    end
end

bandwidths(R::QROperatorR) = (0, bandwidth(R.QR.R_cache,1) + bandwidth(R.QR.R_cache,2))

struct QROperatorQ{QRT,T} <: Operator{T}
    QR::QRT
end


QROperatorQ(QR) = QROperatorQ{typeof(QR),eltype(QR)}(QR)


domainspace(Q::QROperatorQ) = ℓ⁰
rangespace(Q::QROperatorQ) = rangespace(Q.QR)

getindex(Q::QROperatorQ, k::Integer, j::Integer) = (mul_coefficients(Q',eltype(Q)[zeros(k-1);1]))[j]


function getproperty(F::QROperator, d::Symbol)
    if d == :R
        return QROperatorR(F)
    elseif d == :Q
        return QROperatorQ(F)
    else
        getfield(F, d)
    end
end

# iteration for destructuring into components
Base.iterate(S::QROperator) = (S.Q, Val(:R))
Base.iterate(S::QROperator, ::Val{:R}) = (S.R, Val(:done))
Base.iterate(S::QROperator, ::Val{:done}) = nothing


# override for custom data types
QROperator(R::CachedOperator{T,AM}) where {T,AM<:AbstractMatrix} =
    error("Cannot create a QR factorization for $(typeof(R))")


adjoint(Q::QROperatorQ) = Adjoint(Q)
size(Q::Adjoint{<:Any,<:QROperatorQ}) = (size(parent(Q),2), size(parent(Q),1))

function qr!(A::CachedOperator; cached::Int=0)
    QR = QROperator(A)
    if cached ≠ 0
        resizedata!(QR,:,cached)
    end
    QR
end

"""
    qr(A::Operator)

returns a cached QR factorization of the Operator `A`.  The result `QR`
enables solving of linear equations: if `u=QR\b`, then `u`
approximately satisfies `A*u = b`.
"""
function qr(A::Operator; cached::Int=0)
    if isambiguous(domainspace(A)) || isambiguous(rangespace(A))
        throw(ArgumentError("Only non-ambiguous operators can be factorized."))
    end
    qr!(cache(A;padding=true);cached=cached)
end


factorize(A::Operator) = qr(A)

for OP in (:qr, :factorize)
    @eval begin
        $OP(A::AbstractVector{<:Operator}) = $OP(interlace(A))
        $OP(A::AbstractMatrix{<:Operator}) = $OP(interlace(A))
        $OP(A::AbstractArray{<:Operator}) = $OP(interlace(A))
    end
end


function det(R::QROperatorR;maxiterations::Int=10_000)
    QR = R.QR
    RD = R.QR.R_cache
    resizedata!(QR,:,1)
    ret = -RD[1,1]
    k = 2
    while abs(abs(RD[k-1,k-1])-1) > eps(eltype(R))
        resizedata!(QR,:,k)
        ret *= -RD[k,k]
        k+=1
        k > maxiterations && error("Determinant unlikely to converge after 10_000 iterations.")
    end

    ret
end
det(QR::QROperatorQ) = 1

det(QR::QROperator) = det(QR.R)
det(A::Operator) = det(qr(A))




## Multiplication routines


# Q

mul_coefficients(At::Transpose{T,<:QROperatorQ{T}},B::AbstractVector{T}) where {T<:Real} = parent(At)'*B
mul_coefficients(At::Transpose{T,<:QROperatorQ{T}},B::AbstractMatrix{T}) where {T<:Real} = parent(At)'*B

mul_coefficients(Ac::Adjoint{T,<:QROperatorQ{QR,T}},B::AbstractVector{T};tolerance=eps(eltype(Ac))/10,maxlength=1000000) where {QR,T} =
        mulpars(Ac,B,tolerance,maxlength)

mul_coefficients(Ac::Adjoint{T,<:QROperatorQ{QR,T}},B::AbstractVector{V};opts...) where {QR,T,V} =
    mul_coefficients(Ac,AbstractVector{T}(B); opts...)


function *(Ac::Adjoint{<:Any,<:QROperatorQ}, b::AbstractVector; kwds...)
    A = parent(Ac)
    Fun(domainspace(A),mul_coefficients(Ac,coefficients(b,rangespace(A));kwds...))
end

function *(Ac::Adjoint{<:Any,<:QROperatorQ}, b; kwds...)
    A = parent(Ac)
    Fun(domainspace(A),mul_coefficients(Ac,coefficients(b,rangespace(A));kwds...))
end


ldiv_coefficients(A::QROperatorQ, B; opts...) = mul_coefficients(A', B; opts...)
\(A::QROperatorQ, B::Fun; opts...) = *(A', B; opts...)


# R
function ldiv_coefficients(R::QROperatorR, b::AbstractVector)
    if length(b) > R.QR.ncols
        # upper triangularize columns
        resizedata!(R.QR, :, length(b))
    end
    UpperTriangular(view(R.QR.R_cache.data, 1:length(b), 1:length(b))) \ b
end

\(R::QROperatorR,b::Fun{SequenceSpace};kwds...) =
    Fun(domainspace(R),ldiv_coefficients(R,b.coefficients;kwds...))
\(A::QROperatorR,b::Fun;kwds...) = error("\\ not implement for $(typeof(b)) right-hand sides")


# QR

for TYP in (:Real,:Complex,:Number)
    @eval ldiv_coefficients(QR::QROperator{CO,MT,T},b::AbstractVector{T}; kwds...) where {CO,MT,T<:$TYP} =
        ldiv_coefficients(QR.R, mul_coefficients(QR.Q',b;kwds...))
end


function ldiv_coefficients(QR::QROperator{CO,MT,T},b::AbstractVector{V};kwds...) where {CO,MT,T,V<:Number}
    TV = promote_type(T,V)
    ldiv_coefficients(convert(Operator{TV}, QR),convert(Vector{TV}, b);kwds...)
end

function ldiv_coefficients(QR::QROperator{CO,MT,T},b::AbstractVector{V};kwds...) where {CO,MT,T<:Real,V<:Complex}
    a=ldiv_coefficients(QR,real(b);kwds...)
    b=im*ldiv_coefficients(QR,imag(b);kwds...)
    n=max(length(a),length(b))
    pad!(a,n)+pad!(b,n)
end
ldiv_coefficients(QR::QROperator{CO,MT,T},b::AbstractVector{V};kwds...) where {CO,MT,T<:Complex,V<:Real} =
    ldiv_coefficients(QR,Vector{T}(b);kwds...)


\(A::QROperator,b::Fun;kwds...) =
    Fun(domainspace(A),ldiv_coefficients(A,coefficients(b,rangespace(A));kwds...))

\(A::QROperator,B::MatrixFun;kwds...) = \(A,Array(B);kwds...)
