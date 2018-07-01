
recA(::Type{T},::Chebyshev,k) where {T} = 2one(T)
recB(::Type{T},::Chebyshev,_) where {T} = zero(T)
recC(::Type{T},::Chebyshev,k) where {T} = one(T)   # one(T) ensures we get correct type

recα(::Type{T},::Chebyshev,_) where {T} = zero(T)
recβ(::Type{T},::Chebyshev,k) where {T} = ifelse(k==1,one(T),one(T)/2)   # one(T) ensures we get correct type,ifelse ensures inlining
recγ(::Type{T},::Chebyshev,k) where {T} = one(T)/2   # one(T) ensures we get correct type





## Evaluation

Evaluation(S::Chebyshev,x::Bool,o::Integer) =
    ConcreteEvaluation(S,x,o)

Evaluation(S::Chebyshev,x::Number,o::Integer) =
    o==0 ? ConcreteEvaluation(S,x,o) : EvaluationWrapper(S,x,o,Evaluation(x)*Derivative(S,o))

function evaluatechebyshev(n::Integer,x::T) where T<:Number
    if n == 1
        [one(T)]
    else
        p = zeros(T,n)
        p[1] = one(T)
        p[2] = x

        for j=2:n-1
            p[j+1] = 2x*p[j] - p[j-1]
        end

        p
    end
end


function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(first)},j::Integer) where {DD<:Segment,RR}
    T=eltype(op)
    if op.order == 0
        ifelse(isodd(j),  # right rule
            one(T),
            -one(T))
    else
        #TODO: Fast version
        op[j:j][1]
    end
end

function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(last)},j::Integer) where {DD<:Segment,RR}
    T=eltype(op)
    if op.order == 0
        one(T)
    else
        #TODO: Fast version
        op[j:j][1]
    end
end



function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(first)},k::AbstractRange) where {DD<:Segment,RR}
    T=eltype(op)
    x = op.x
    d = domain(op)
    p = op.order
    cst = T((2/(d.b-d.a))^p)
    n=length(k)

    ret = Array{T}(undef, n)
    k1=1-first(k)
    @simd for j=k
        @inbounds ret[j+k1]=(-1)^(p+1)*(-one(T))^j
    end

    for m=0:p-1
        k1=1-first(k)
        @simd for j=k
            @inbounds ret[j+k1] *= (j-1)^2-m^2
        end
        scal!(T(1/(2m+1)), ret)
    end

    scal!(cst,ret)
end

function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(last)},k::AbstractRange) where {DD<:Segment,RR}
    T=eltype(op)
    x = op.x
    d = domain(op)
    p = op.order
    cst = T((2/(d.b-d.a))^p)
    n=length(k)

    ret = fill(one(T),n)

    for m=0:p-1
        k1=1-first(k)
        @simd for j=k
            @inbounds ret[j+k1] *= (j-1)^2-m^2
        end
        scal!(T(1/(2m+1)), ret)
    end

    scal!(cst,ret)
end

function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},M,OT,T},
                j::Integer) where {DD<:Segment,RR,M<:Real,OT,T}
    if op.order == 0
        T(evaluatechebyshev(j,tocanonical(domain(op),op.x))[end])
    else
        error("Only zero–second order implemented")
    end
end

function getindex(op::ConcreteEvaluation{Chebyshev{DD,RR},M,OT,T},
                                             k::AbstractRange) where {DD<:Segment,RR,M<:Real,OT,T}
    if op.order == 0
        Array{T}(evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k])
    else
        error("Only zero–second order implemented")
    end
end

function Dirichlet(S::Chebyshev,order)
    order == 0 && return ConcreteDirichlet(S,ArraySpace([ConstantSpace(Point(first(domain(S)))),
                                                         ConstantSpace(Point(last(domain(S))))]),0)
    default_Dirichlet(S,order)
end


function getindex(op::ConcreteDirichlet{<:Chebyshev},
                                             k::Integer,j::Integer)
    if op.order == 0
        k == 1 && iseven(j) && return -one(eltype(op))
        return one(eltype(op))
    else
        error("Only zero Dirichlet conditions implemented")
    end
end

function convert(::Type{Matrix},
                 S::SubOperator{T,ConcreteDirichlet{C,V,T},
                                Tuple{UnitRange{Int},UnitRange{Int}}}) where {C<:Chebyshev,V,T}
    ret = Array{T}(undef, size(S)...)
    kr,jr = parentindices(S)
    isempty(kr) && return ret
    isempty(jr) && return ret
    if first(kr) == 1
        if isodd(jr[1])
            ret[1,1:2:end] = one(T)
            ret[1,2:2:end] = -one(T)
        else
            ret[1,1:2:end] = -one(T)
            ret[1,2:2:end] = one(T)
        end
    end
    if last(kr) == 2
        ret[end,:] = one(T)
    end
    return ret
end
#

# Multiplication

Base.stride(M::ConcreteMultiplication{U,V}) where {U<:Chebyshev,V<:Chebyshev} =
    stride(M.f)


function chebmult_getindex(cfs::AbstractVector,k::Integer,j::Integer)
    n=length(cfs)

    ret = zero(eltype(cfs))

    n == 0 && return ret

    # Toeplitz part
    if k == j
        ret += cfs[1]
    elseif k > j && k-j+1 ≤ n
        ret += cfs[k-j+1]/2
    elseif k < j && j-k+1 ≤ n
        ret += cfs[j-k+1]/2
    end

    # Hankel part
    if k ≥ 2 && k+j-1 ≤ n
        ret += cfs[k+j-1]/2
    end

    ret
end


getindex(M::ConcreteMultiplication{C,C,T},k::Integer,j::Integer) where {T,C<:Chebyshev} =
    chebmult_getindex(coefficients(M.f),k,j)

getindex(M::ConcreteMultiplication{C,PS,T},k::Integer,j::Integer) where {PS<:PolynomialSpace,T,C<:Chebyshev} =
    M[k:k,j:j][1,1]


function convert(::Type{BandedMatrix},S::SubOperator{T,ConcreteMultiplication{C,C,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {C<:Chebyshev,T}
    ret = BandedMatrix(Zeros, S)

    kr,jr=parentindices(S)
    cfs=parent(S).f.coefficients

    isempty(cfs) && return ret

    # Toeplitz part
    sym_toeplitz_axpy!(1.0,0.5,cfs,kr,jr,ret)

    #Hankel part
    hankel_axpy!(0.5,cfs,kr,jr,ret)

    # divide first row by half
    if first(kr)==1
        if first(jr)==1
            ret[1,1]+=0.5cfs[1]
        end

        for j=1:min(1+ret.u,size(ret,2))
            ret[1,j]/=2
        end
    end


    ret
end



## Derivative

Derivative(sp::Chebyshev{DD},order::Integer) where {DD<:Segment} =
    ConcreteDerivative(sp,order)


rangespace(D::ConcreteDerivative{Chebyshev{DD,RR}}) where {DD<:Segment,RR} =
    Ultraspherical(D.order,domain(D))
bandinds(D::ConcreteDerivative{Chebyshev{DD,RR}}) where {DD<:Segment,RR} = D.order,D.order
Base.stride(D::ConcreteDerivative{Chebyshev{DD,RR}}) where {DD<:Segment,RR} = D.order

function getindex(D::ConcreteDerivative{Chebyshev{DD,RR},K,T},k::Integer,j::Integer) where {DD<:Segment,RR,K,T}
    m=D.order
    d=domain(D)

    if j==k+m
        C=T(pochhammer(one(T),m-1)/2*(4/(d.b-d.a))^m)
        T(C*(m+k-one(T)))
    else
        zero(T)
    end
end

linesum(f::Fun{Chebyshev{DD,RR}}) where {DD<:Segment,RR} =
    sum(setcanonicaldomain(f))*arclength(domain(f))/2



## Clenshaw-Curtis functional

for (Func,Len) in ((:DefiniteIntegral,:complexlength),(:DefiniteLineIntegral,:arclength))
    ConcFunc = Meta.parse("Concrete"*string(Func))
    @eval begin
        $Func(S::Chebyshev{D}) where {D<:Segment} = $ConcFunc(S)
        function getindex(Σ::$ConcFunc{Chebyshev{D,R},T},k::Integer) where {D<:Segment,R,T}
            d = domain(Σ)
            C = $Len(d)/2

            isodd(k) ? T(2C/(k*(2-k))) : zero(T)
        end
    end
end



ReverseOrientation(S::Chebyshev) = ReverseOrientationWrapper(NegateEven(S,reverseorientation(S)))
Reverse(S::Chebyshev) = ReverseWrapper(NegateEven(S,S))
