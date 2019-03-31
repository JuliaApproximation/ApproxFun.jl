
export Chebyshev, NormalizedChebyshev

"""
`Chebyshev()` is the space spanned by the Chebyshev polynomials
```
    T_0(x),T_1(x),T_2(x),…
```
where `T_k(x) = cos(k*acos(x))`.  This is the default space
as there exists a fast transform and general smooth functions on `[-1,1]`
can be easily resolved.
"""
struct Chebyshev{D<:Domain,R} <: PolynomialSpace{D,R}
    domain::D
    function Chebyshev{D,R}(d) where {D,R}
        isempty(d) && throw(ArgumentError("Domain cannot be empty"))
        new(d)
    end
    Chebyshev{D,R}() where {D,R} = new(convert(D, ChebyshevInterval()))
end

Chebyshev(d::Domain) = Chebyshev{typeof(d),real(prectype(d))}(d)
Chebyshev() = Chebyshev(ChebyshevInterval())
Chebyshev(d) = Chebyshev(Domain(d))
NormalizedChebyshev() = NormalizedPolynomialSpace(Chebyshev())
NormalizedChebyshev(d) = NormalizedPolynomialSpace(Chebyshev(d))

normalization(::Type{T}, sp::Chebyshev, k::Int) where T = T(π)/(2-FastTransforms.δ(k,0))

Space(d::SegmentDomain) = Chebyshev(d)
function Space(d::AbstractInterval)
    a,b = endpoints(d)
    if isinf(norm(a)) && isinf(norm(b))
        Chebyshev(Line(d))
    elseif isinf(norm(a)) || isinf(norm(b))
        Chebyshev(Ray(d))
    else
        Chebyshev(d)
    end
end


setdomain(S::Chebyshev, d::Domain) = Chebyshev(d)

ones(::Type{T}, S::Chebyshev) where {T<:Number} = Fun(S,fill(one(T),1))
ones(S::Chebyshev) = Fun(S,fill(1.0,1))

function Base.first(f::Fun{<:Chebyshev})
    n = ncoefficients(f)
    n == 0 && return zero(cfstype(f))
    n == 1 && return f.coefficients[1]
    foldr(-,coefficients(f))
end

Base.last(f::Fun{<:Chebyshev}) = reduce(+,coefficients(f))

spacescompatible(a::Chebyshev,b::Chebyshev) = domainscompatible(a,b)
hasfasttransform(::Chebyshev) = true


## Transform

transform(::Chebyshev,vals::AbstractVector,plan) = plan*vals
itransform(::Chebyshev,cfs::AbstractVector,plan) = plan*cfs
plan_transform(::Chebyshev,vals::AbstractVector) = plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::AbstractVector) = plan_ichebyshevtransform(cfs)

## Evaluation

clenshaw(sp::Chebyshev, c::AbstractVector, x::AbstractArray) =
    clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))

clenshaw(::Chebyshev,c::AbstractVector,x) = chebyshev_clenshaw(c,x)

clenshaw(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S<:Chebyshev,V}=
    clenshaw(c,collect(x),plan)

#TODO: This modifies x, which is not threadsafe
clenshaw(c::AbstractVector,x::Vector,plan::ClenshawPlan{<:Chebyshev}) = chebyshev_clenshaw(c,x,plan)

function clenshaw(c::AbstractMatrix{T},x::T,plan::ClenshawPlan{S,T}) where {S<:Chebyshev,T<:Number}
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    m,n=size(c) # m is # of coefficients, n is # of funs

    @inbounds for i = 1:n
        bk1[i] = zero(T)
        bk2[i] = zero(T)
    end
    x = 2x

    @inbounds for k=m:-1:2
        for j=1:n
            ck = c[k,j]
            bk[j] = muladd(x,bk1[j],ck - bk2[j])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    x = x/2
    @inbounds for i = 1:n
        ce = c[1,i]
        bk[i] = muladd(x,bk1[i],ce - bk2[i])
    end

    bk
end

function clenshaw(c::AbstractMatrix{T},x::AbstractVector{T},plan::ClenshawPlan{S,T}) where {S<:Chebyshev,T<:Number}
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    m,n=size(c) # m is # of coefficients, n is # of funs

    @inbounds for i = 1:n
        x[i] = 2x[i]
        bk1[i] = zero(T)
        bk2[i] = zero(T)
    end

    @inbounds for k=m:-1:2
        for j=1:n
            ck = c[k,j]
            bk[j] = muladd(x[j],bk1[j],ck - bk2[j])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end


    @inbounds for i = 1:n
        x[i] = x[i]/2
        ce = c[1,i]
        bk[i] = muladd(x[i],bk1[i],ce - bk2[i])
    end

    bk
end

# overwrite x

function clenshaw!(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S<:Chebyshev,V}
    N,n = length(c),length(x)

    if isempty(c)
        fill!(x,0)
        return x
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    @inbounds for i = 1:n
        x[i] = 2x[i]
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck = c[k]
        for i = 1:n
            bk[i] = muladd(x[i],bk1[i],ck - bk2[i])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ce = c[1]
    @inbounds for i = 1:n
        x[i] = x[i]/2
        x[i] = muladd(x[i],bk1[i],ce-bk2[i])
    end

    x
end

## Calculus


# diff T -> U, then convert U -> T
integrate(f::Fun{Chebyshev{D,R}}) where {D<:IntervalOrSegment,R} =
    Fun(f.space,fromcanonicalD(f,0)*ultraint!(ultraconversion(f.coefficients)))
differentiate(f::Fun{Chebyshev{D,R}}) where {D<:IntervalOrSegment,R} =
    Fun(f.space,1/fromcanonicalD(f,0)*ultraiconversion(ultradiff(f.coefficients)))





## Multivariate

# determine correct parameter to have at least
# N point
_padua_length(N) = Int(cld(-3+sqrt(1+8N),2))

function squarepoints(::Type{T}, N) where T
    pts = paduapoints(T, _padua_length(N))
    n = size(pts,1)
    ret = Array{Vec{2,T}}(undef, n)
    @inbounds for k=1:n
        ret[k] = Vec{2,T}(pts[k,1],pts[k,2])
    end
    ret
end

points(S::TensorSpace{<:Tuple{<:Chebyshev{<:ChebyshevInterval},<:Chebyshev{<:ChebyshevInterval}}}, N) =
    squarepoints(real(prectype(S)), N)

function points(S::TensorSpace{<:Tuple{<:Chebyshev,<:Chebyshev}},N)
    T = real(prectype(S))
    pts = squarepoints(T, N)
    pts .= fromcanonical.(Ref(domain(S)), pts)
    pts
end

plan_transform(S::TensorSpace{<:Tuple{<:Chebyshev,<:Chebyshev}},v::AbstractVector) =
    plan_paduatransform!(v,Val{false})

transform(S::TensorSpace{<:Tuple{<:Chebyshev,<:Chebyshev}},v::AbstractVector,
        plan=plan_transform(S,v)) = plan*copy(v)

plan_itransform(S::TensorSpace{<:Tuple{<:Chebyshev,<:Chebyshev}},v::AbstractVector) =
     plan_ipaduatransform!(eltype(v),sum(1:nblocks(Fun(S,v))),Val{false})

itransform(S::TensorSpace{<:Tuple{<:Chebyshev,<:Chebyshev}},v::AbstractVector,
         plan=plan_itransform(S,v)) = plan*pad(v,sum(1:nblocks(Fun(S,v))))


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::ProductFun{<:Chebyshev,<:Chebyshev}) =
        ProductFun(chebyshevtransform($op.(values(f))),space(f))
end



reverseorientation(f::Fun{<:Chebyshev}) =
    Fun(Chebyshev(reverseorientation(domain(f))),alternatesign!(copy(f.coefficients)))


include("ChebyshevOperators.jl")
