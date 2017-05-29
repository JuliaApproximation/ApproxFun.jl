
recA{T}(::Type{T},::Chebyshev,k) = 2one(T)
recB{T}(::Type{T},::Chebyshev,::) = zero(T)
recC{T}(::Type{T},::Chebyshev,k) = one(T)   # one(T) ensures we get correct type

recα{T}(::Type{T},::Chebyshev,::) = zero(T)
recβ{T}(::Type{T},::Chebyshev,k) = ifelse(k==1,one(T),one(T)/2)   # one(T) ensures we get correct type,ifelse ensures inlining
recγ{T}(::Type{T},::Chebyshev,k) = one(T)/2   # one(T) ensures we get correct type





## Evaluation

Evaluation(S::Chebyshev,x::Bool,o::Integer) =
    ConcreteEvaluation(S,x,o)

Evaluation(S::Chebyshev,x::Number,o::Integer) =
    o==0?ConcreteEvaluation(S,x,o):EvaluationWrapper(S,x,o,Evaluation(x)*Derivative(S,o))

function evaluatechebyshev{T<:Number}(n::Integer,x::T)
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


function getindex{DD<:Segment,RR}(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(first)},j::Integer)
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

function getindex{DD<:Segment,RR}(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(last)},j::Integer)
    T=eltype(op)
    if op.order == 0
        one(T)
    else
        #TODO: Fast version
        op[j:j][1]
    end
end



function getindex{DD<:Segment,RR}(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(first)},k::Range)
    T=eltype(op)
    x = op.x
    d = domain(op)
    p = op.order
    cst = T((2/(d.b-d.a))^p)
    n=length(k)

    ret = Array{T}(n)
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

function getindex{DD<:Segment,RR}(op::ConcreteEvaluation{Chebyshev{DD,RR},typeof(last)},k::Range)
    T=eltype(op)
    x = op.x
    d = domain(op)
    p = op.order
    cst = T((2/(d.b-d.a))^p)
    n=length(k)

    ret = ones(T,n)

    for m=0:p-1
        k1=1-first(k)
        @simd for j=k
            @inbounds ret[j+k1] *= (j-1)^2-m^2
        end
        scal!(T(1/(2m+1)), ret)
    end

    scal!(cst,ret)
end

function getindex{DD<:Segment,RR,M<:Real,OT,T}(op::ConcreteEvaluation{Chebyshev{DD,RR},M,OT,T},
                                             j::Integer)
    if op.order == 0
        T(evaluatechebyshev(j,tocanonical(domain(op),op.x))[end])
    else
        error("Only zero–second order implemented")
    end
end

function getindex{DD<:Segment,RR,M<:Real,OT,T}(op::ConcreteEvaluation{Chebyshev{DD,RR},M,OT,T},
                                             k::Range)
    if op.order == 0
        Array{T}(evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k])
    else
        error("Only zero–second order implemented")
    end
end


# Multiplication

Base.stride{U<:Chebyshev,V<:Chebyshev}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)


function chebmult_getindex(cfs::AbstractVector,k::Integer,j::Integer)
    n=length(cfs)

    ret=zero(eltype(cfs))

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


getindex{T,C<:Chebyshev}(M::ConcreteMultiplication{C,C,T},k::Integer,j::Integer) =
    chebmult_getindex(coefficients(M.f),k,j)

getindex{PS<:PolynomialSpace,T,C<:Chebyshev}(M::ConcreteMultiplication{C,PS,T},k::Integer,j::Integer) =
    M[k:k,j:j][1,1]


function Base.convert{C<:Chebyshev,T}(::Type{BandedMatrix},S::SubOperator{T,ConcreteMultiplication{C,C,T},Tuple{UnitRange{Int},UnitRange{Int}}})
    ret=bzeros(S)

    kr,jr=parentindexes(S)
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

Derivative{DD<:Segment}(sp::Chebyshev{DD},order::Integer) =
    ConcreteDerivative(sp,order)


rangespace{DD<:Segment,RR}(D::ConcreteDerivative{Chebyshev{DD,RR}}) =
    Ultraspherical(D.order,domain(D))
bandinds{DD<:Segment,RR}(D::ConcreteDerivative{Chebyshev{DD,RR}}) = D.order,D.order
Base.stride{DD<:Segment,RR}(D::ConcreteDerivative{Chebyshev{DD,RR}}) = D.order

function getindex{DD<:Segment,RR,K,T}(D::ConcreteDerivative{Chebyshev{DD,RR},K,T},k::Integer,j::Integer)
    m=D.order
    d=domain(D)

    if j==k+m
        C=T(pochhammer(one(T),m-1)/2*(4/(d.b-d.a))^m)
        T(C*(m+k-one(T)))
    else
        zero(T)
    end
end

linesum{DD<:Segment,RR}(f::Fun{Chebyshev{DD,RR}}) =
    sum(setcanonicaldomain(f))*arclength(domain(f))/2



## Clenshaw-Curtis functional

for (Func,Len) in ((:DefiniteIntegral,:complexlength),(:DefiniteLineIntegral,:arclength))
    ConcFunc = parse("Concrete"*string(Func))
    @eval begin
        $Func{D<:Segment}(S::Chebyshev{D}) = $ConcFunc(S)
        function getindex{D<:Segment,R,T}(Σ::$ConcFunc{Chebyshev{D,R},T},k::Integer)
            d = domain(Σ)
            C = $Len(d)/2

            isodd(k) ? T(2C/(k*(2-k))) : zero(T)
        end
    end
end



ReverseOrientation(S::Chebyshev) = ReverseOrientationWrapper(SpaceOperator(NegateEven(),S,reverseorientation(S)))
Reverse(S::Chebyshev) = ReverseWrapper(SpaceOperator(NegateEven(),S,S))
