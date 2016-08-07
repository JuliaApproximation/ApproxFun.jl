

##
# Represent an m x n rank-r matrix
# A = U V^T
##

export LowRankMatrix

type LowRankMatrix{T} <: AbstractMatrix{T}
    U::Matrix{T} # m x r Matrix
    V::Matrix{T} # n x r Matrix

    function LowRankMatrix(U::Matrix{T},V::Matrix{T})
        m,r = size(U)
        n,rv = size(V)
        @assert r == rv
        new(U,V)
    end
end

LowRankMatrix{T}(U::Matrix{T},V::Matrix{T})=LowRankMatrix{T}(U,V)

LowRankMatrix(U::Matrix,V::Matrix) = LowRankMatrix{promote_type(eltype(U),eltype(V))}(promote(U,V)...)
LowRankMatrix(U::Vector,V::Matrix) = LowRankMatrix(reshape(U,length(U),1),V)
LowRankMatrix(U::Matrix,V::Vector) = LowRankMatrix(U,reshape(V,length(V),1))
LowRankMatrix(U::Vector,V::Vector) = LowRankMatrix(reshape(U,length(U),1),reshape(V,length(V),1))
LowRankMatrix(a::Number,m::Int,n::Int) = LowRankMatrix(a*ones(eltype(a),m),ones(eltype(a),n))

LowRankMatrix{T}(::Type{T},m::Int,n::Int,r::Int) = LowRankMatrix(Array(T,m,r),Array(T,n,r))
lrzeros{T}(::Type{T},m::Int,n::Int,r::Int) = LowRankMatrix(zeros(T,m,r),zeros(T,n,r))

Base.similar(L::LowRankMatrix, T, dims::Dims) = (@assert length(dims) == 2;r = rank(L); LowRankMatrix(Array(T,dims[1],r),Array(T,dims[2],r)))
Base.similar{T}(L::LowRankMatrix{T}) = ((m,n) = size(L); r = rank(L); LowRankMatrix(Array(T,m,r),Array(T,n,r)))
Base.similar{T}(L::LowRankMatrix{T}, dims::Dims) = (@assert length(dims) == 2;r = rank(L); LowRankMatrix(Array(T,dims[1],r),Array(T,dims[2],r)))
Base.similar{T}(L::LowRankMatrix{T}, m::Int) = Array(T, m)
Base.similar{T}(L::LowRankMatrix{T}, S) = ((m,n) = size(L); r = rank(L); LowRankMatrix(Array(S,m,r),Array(S,n,r)))

function LowRankMatrix(A::Matrix)
    U,Σ,V = svd(A)
    r = refactorsvd!(U,Σ,V)
    LowRankMatrix(U[:,1:r],V[:,1:r])
end

balance!{T<:Union{Integer,Rational}}(U::Matrix{T},V::Matrix{T},m::Int,n::Int,r::Int) = U,V
function balance!{T}(U::Matrix{T},V::Matrix{T},m::Int,n::Int,r::Int)
    for k=1:r
        uk = zero(T)
        for i=1:m
            @inbounds uk += abs2(U[i,k])
        end
        vk = zero(T)
        for j=1:n
            @inbounds vk += abs2(V[j,k])
        end
        uk,vk = sqrt(uk),sqrt(vk)
        σk = sqrt(uk*vk)
        if abs2(uk) ≥ eps(T)^2 && abs2(vk) ≥ eps(T)^2
            uk,vk = σk/uk,σk/vk
            for i=1:m
                @inbounds U[i,k] *= uk
            end
            for j=1:n
                @inbounds V[j,k] *= vk
            end
        end
    end
    U,V
end

function refactorsvd!{S,T}(U::Matrix{S},Σ::Vector{T},V::Matrix{S})
    conj!(V)
    σmax = Σ[1]
    r=max(1,count(s->s>10σmax*eps(T),Σ))
    m,n = size(U,1),size(V,1)
    for k=1:r
        σk = sqrt(Σ[k])
        for i=1:m
            @inbounds U[i,k] *= σk
        end
        for j=1:n
            @inbounds V[j,k] *= σk
        end
    end
    r
end

Base.convert{T}(::Type{LowRankMatrix{T}},L::LowRankMatrix) = LowRankMatrix{T}(convert(Matrix{T},L.U),convert(Matrix{T},L.V))
Base.convert{T}(::Type{Matrix{T}},L::LowRankMatrix) = convert(Matrix{T},full(L))
Base.promote_rule{T,V}(::Type{LowRankMatrix{T}},::Type{LowRankMatrix{V}})=LowRankMatrix{promote_type(T,V)}
Base.promote_rule{T,V}(::Type{LowRankMatrix{T}},::Type{Matrix{V}})=Matrix{promote_type(T,V)}

Base.size(L::LowRankMatrix) = size(L.U,1),size(L.V,1)
Base.rank(L::LowRankMatrix) = size(L.U,2)
Base.transpose(L::LowRankMatrix) = LowRankMatrix(L.V,L.U)
Base.ctranspose{T<:Real}(L::LowRankMatrix{T}) = LowRankMatrix(L.V,L.U)
Base.ctranspose(L::LowRankMatrix) = LowRankMatrix(conj(L.V),conj(L.U))
Base.fill!{T}(L::LowRankMatrix{T}, x::T) = (fill!(L.U, sqrt(abs(x)/rank(L)));fill!(L.V,sqrt(abs(x)/rank(L))/sign(x)); L)

function Base.getindex(L::LowRankMatrix,i::Int,j::Int)
    m,n = size(L)
    if 1 ≤ i ≤ m && 1 ≤ j ≤ n
        ret = zero(eltype(L))
        for k=1:rank(L)
            @inbounds ret = muladd(L.U[i,k],L.V[j,k],ret)
        end
        return ret
    else
        throw(BoundsError())
    end
end
Base.getindex(L::LowRankMatrix,i::Int,jr::Range) = eltype(L)[L[i,j] for j=jr].'
Base.getindex(L::LowRankMatrix,ir::Range,j::Int) = eltype(L)[L[i,j] for i=ir]
Base.getindex(L::LowRankMatrix,ir::Range,jr::Range) = eltype(L)[L[i,j] for i=ir,j=jr]
Base.full(L::LowRankMatrix)=L[1:size(L,1),1:size(L,2)]

# constructors

for op in (:zeros,:eye,:ones,:rand)
    lrop = parse("lr"*string(op))
    @eval begin
        export $lrop
        $lrop(T::Type,m::Int,n::Int) = LowRankMatrix($op(T,m,n))
        $lrop(T::Type,n::Int) = $lrop(T,n,n)
        $lrop(m::Int,n::Int) = LowRankMatrix($op(Float64,m,n))
        $lrop(n::Int) = $lrop(n,n)
    end
end
export lrrandn
lrrandn(::Type{Float64},m::Int,n::Int) = LowRankMatrix(randn(m,n))
lrrandn(::Type{Float64},n::Int) = lrrandn(n,n)
lrrandn(m::Int,n::Int) = lrrandn(m,n)
lrrandn(n::Int) = lrrandn(n,n)

Base.copy(L::LowRankMatrix) = LowRankMatrix(copy(L.U),copy(L.V))
Base.copy!(L::LowRankMatrix,N::LowRankMatrix) = (copy!(L.U,N.U);copy!(L.V,N.V);L)


function pad!(L::LowRankMatrix,n::Integer,::Colon)
    L.U=pad(L.U,n,:)
    L
end
function pad!(L::LowRankMatrix,::Colon,m::Integer)
    L.V=pad(L.V,m,:)
    L
end
pad!(L::LowRankMatrix,n::Integer,m::Integer) = pad!(pad!(L,n,:),:,m)

# algebra

for op in (:+,:-,:.+,:.-)
    @eval begin
        $op(L::LowRankMatrix) = LowRankMatrix($op(L.U),L.V)

        $op(a::Bool,L::LowRankMatrix{Bool}) = error("Not callable")
        $op(L::LowRankMatrix{Bool},a::Bool) = error("Not callable")
        $op(a::Number,L::LowRankMatrix) = $op(LowRankMatrix(a,size(L)...),L)
        $op(L::LowRankMatrix,a::Number) = $op(L,LowRankMatrix(a,size(L)...))

        function $op(L::LowRankMatrix,M::LowRankMatrix)
            @assert size(L) == size(M)
            LowRankMatrix(hcat(L.U,$op(M.U)),hcat(L.V,M.V))
        end
        $op(L::LowRankMatrix,A::Matrix) = $op(promote(L,A)...)
        $op(A::Matrix,L::LowRankMatrix) = $op(promote(A,L)...)
    end
end

*(a::Number,L::LowRankMatrix) = LowRankMatrix(a*L.U,L.V)
*(L::LowRankMatrix,a::Number) = LowRankMatrix(L.U,L.V*a)
.*(a::Number,L::LowRankMatrix) = a*L
.*(L::LowRankMatrix,a::Number) = L*a

function Base.A_mul_B!(b::AbstractVector,L::LowRankMatrix,x::AbstractVector)
    temp = zeros(promote_type(eltype(L),eltype(x)),rank(L))
    At_mul_B!(temp,L.V,x)
    A_mul_B!(b,L.U,temp)
    b
end

function *(L::LowRankMatrix,M::LowRankMatrix)
    T = promote_type(eltype(L),eltype(M))
    temp = zeros(T,rank(L),rank(M))
    At_mul_B!(temp,L.V,M.U)
    V = zeros(T,size(M,2),rank(L))
    A_mul_Bt!(V,M.V,temp)
    LowRankMatrix(copy(L.U),V)
end

function *(L::LowRankMatrix,A::AbstractMatrix)
    V = zeros(promote_type(eltype(L),eltype(A)),size(A,2),rank(L))
    At_mul_B!(V,A,L.V)
    LowRankMatrix(copy(L.U),V)
end

function *(A::AbstractMatrix,L::LowRankMatrix)
    U = zeros(promote_type(eltype(A),eltype(L)),size(A,1),rank(L))
    At_mul_B!(U,A,L.U)
    LowRankMatrix(U,copy(L.V))
end

\(L::LowRankMatrix,b::AbstractVecOrMat) = full(L)\b
