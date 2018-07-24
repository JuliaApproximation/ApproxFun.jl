
export Chebyshev

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
    Chebyshev{D,R}(d) where {D,R} = new(d)
    Chebyshev{D,R}() where {D,R} = new(Segment())
end

Chebyshev(d::Domain) = Chebyshev{typeof(d),real(prectype(d))}(d)
Chebyshev() = Chebyshev(ChebyshevInterval())
Chebyshev(d) = Chebyshev(Domain(d))


Space(d::Segment) = Chebyshev(d)
Space(d::Domains.AbstractInterval) = Chebyshev(d)


setdomain(S::Chebyshev,d::Domain) = Chebyshev(d)

ones(::Type{T},S::Chebyshev) where {T<:Number} = Fun(S,fill(one(T),1))
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


function coefficients(g::AbstractVector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function coefficients(g::AbstractVector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::AbstractVector,plan) = plan*vals
itransform(::Chebyshev,cfs::AbstractVector,plan) = plan*cfs
plan_transform(::Chebyshev,vals::AbstractVector) = plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::AbstractVector) = plan_ichebyshevtransform(cfs)

## Evaluation

clenshaw(sp::Chebyshev,c::AbstractVector,x::AbstractArray) =
    clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))

function clenshaw(::Chebyshev,c::AbstractVector,x)
    N,T = length(c),promote_type(eltype(c),typeof(x))
    if N == 0
        return zero(x)
    elseif N == 1 # avoid issues with NaN x
        return first(c)*one(x)
    end

    x = 2x
    bk1,bk2 = zero(T),zero(T)
    @inbounds for k = N:-1:2
        bk2, bk1 = bk1, muladd(x,bk1,c[k]-bk2)
    end

    muladd(x/2,bk1,c[1]-bk2)
end

clenshaw(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S<:Chebyshev,V}=
    clenshaw(c,collect(x),plan)

#TODO: This modifies x, which is not threadsafe
function clenshaw(c::AbstractVector,x::Vector,plan::ClenshawPlan{S,V}) where {S<:Chebyshev,V}
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(V,n)
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
            bk[i] = muladd(x[i],bk1[i],ck-bk2[i])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ck = c[1]
    @inbounds for i = 1:n
        x[i] = x[i]/2
        bk[i] = muladd(x[i],bk1[i],ck-bk2[i])
    end

    bk
end

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

function points(S::TensorSpace{Tuple{Chebyshev{D,R},Chebyshev{D,R}}},N) where {D,R}
    if domain(S) == Segment()^2
        pts=paduapoints(real(prectype(D)),Int(cld(-3+sqrt(1+8N),2)))
        T=eltype(pts)
        ret=Array{Vec{2,T}}(undef, size(pts,1))
        @inbounds for k in eachindex(ret)
            ret[k]=Vec{2,T}(pts[k,1],pts[k,2])
        end
        ret
    else
        fromcanonical.(Ref(S),points(Chebyshev()^2,N))
    end
end

plan_transform(S::TensorSpace{Tuple{Chebyshev{D,R},Chebyshev{D,R}}},v::AbstractVector) where {D,R} =
    plan_paduatransform!(v,Val{false})

transform(S::TensorSpace{Tuple{Chebyshev{D,R},Chebyshev{D,R}}},v::AbstractVector,
        plan=plan_transform(S,v)) where {D,R} = plan*copy(v)

plan_itransform(S::TensorSpace{Tuple{Chebyshev{D,R},Chebyshev{D,R}}},v::AbstractVector) where {D,R} =
     plan_ipaduatransform!(eltype(v),sum(1:nblocks(Fun(S,v))),Val{false})

itransform(S::TensorSpace{Tuple{Chebyshev{D,R},Chebyshev{D,R}}},v::AbstractVector,
         plan=plan_itransform(S,v)) where {D,R} = plan*pad(v,sum(1:nblocks(Fun(S,v))))


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::ProductFun{S,V}) where {S<:Chebyshev,V<:Chebyshev} =
        ProductFun(chebyshevtransform($op.(values(f))),space(f))
end



reverseorientation(f::Fun{C}) where {C<:Chebyshev} =
    Fun(Chebyshev(reverse(domain(f))),alternatesign!(copy(f.coefficients)))


include("ChebyshevOperators.jl")
