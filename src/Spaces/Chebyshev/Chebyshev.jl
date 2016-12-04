
export Chebyshev

doc"""
`Chebyshev()` is the space spanned by the Chebyshev polynomials
```
    T_0(x),T_1(x),T_2(x),…
```
where `T_k(x) = cos(k*acos(x))`.  This is the default space
as there exists a fast transform and general smooth functions on `[-1,1]`
can be easily resolved.
"""
immutable Chebyshev{D<:Domain} <: PolynomialSpace{D}
    domain::D
    Chebyshev(d) = new(d)
    Chebyshev() = new(Segment())
end

Chebyshev() = Chebyshev{Segment{Float64}}()
Chebyshev(d::Domain) = Chebyshev{typeof(d)}(d)
Chebyshev(d) = Chebyshev(Domain(d))


Space(d::Segment) = Chebyshev(d)


setdomain(S::Chebyshev,d::Domain) = Chebyshev(d)

Base.ones{T<:Number}(::Type{T},S::Chebyshev) = Fun(S,ones(T,1))
Base.ones(S::Chebyshev) = Fun(S,ones(1))

Base.first{D}(f::Fun{Chebyshev{D}}) = foldr(-,coefficients(f))
Base.last{D}(f::Fun{Chebyshev{D}}) = reduce(+,coefficients(f))
identity_fun(d::Chebyshev) = identity_fun(domain(d))

spacescompatible(a::Chebyshev,b::Chebyshev) = domainscompatible(a,b)
hasfasttransform(::Chebyshev) = true


function coefficients(g::Vector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function coefficients(g::Vector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::Vector,plan) = plan*vals
itransform(::Chebyshev,cfs::Vector,plan) = plan*cfs
plan_transform(::Chebyshev,vals::Vector) = plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::Vector) = plan_ichebyshevtransform(cfs)

## Evaluation

clenshaw(sp::Chebyshev,c::AbstractVector,x::AbstractArray) =
    clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))

function clenshaw(::Chebyshev,c::AbstractVector,x)
    N,T = length(c),promote_type(eltype(c),typeof(x))
    if isempty(c)
        return zero(x)
    end

    x = 2x
    bk1,bk2 = zero(T),zero(T)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd(x,bk1,c[k]-bk2)
    end

    muladd(x/2,bk1,c[1]-bk2)
end

clenshaw{S<:Chebyshev,V}(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V})=
    clenshaw(c,collect(x),plan)

#TODO: This modifies x, which is not threadsafe
function clenshaw{S<:Chebyshev,V}(c::AbstractVector,x::Vector,plan::ClenshawPlan{S,V})
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

function clenshaw{S<:Chebyshev,T<:Number}(c::AbstractMatrix{T},x::T,plan::ClenshawPlan{S,T})
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

function clenshaw{S<:Chebyshev,T<:Number}(c::AbstractMatrix{T},x::AbstractVector{T},plan::ClenshawPlan{S,T})
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

function clenshaw!{S<:Chebyshev,V}(c::Vector,x::Vector,plan::ClenshawPlan{S,V})
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
integrate{D<:Segment}(f::Fun{Chebyshev{D}}) =
    Fun(f.space,fromcanonicalD(f,0)*ultraint!(ultraconversion(f.coefficients)))
differentiate{D<:Segment}(f::Fun{Chebyshev{D}}) =
    Fun(f.space,1/fromcanonicalD(f,0)*ultraiconversion(ultradiff(f.coefficients)))

## identity_fun





## Multivariate

function points{D}(S::TensorSpace{Tuple{Chebyshev{D},Chebyshev{D}}},N)
    if domain(S) == Segment()^2
        pts=paduapoints(real(eltype(eltype(D))),Int(cld(-3+sqrt(1+8N),2)))
        T=eltype(pts)
        ret=Array(Vec{2,T},size(pts,1))
        @inbounds for k in eachindex(ret)
            ret[k]=Vec{2,T}(pts[k,1],pts[k,2])
        end
        ret
    else
        fromcanonical.(S,points(Chebyshev()^2,N))
    end
end

plan_transform{D}(S::TensorSpace{Tuple{Chebyshev{D},Chebyshev{D}}},v::Vector) =
    plan_paduatransform!(v,Val{false})

transform{D}(S::TensorSpace{Tuple{Chebyshev{D},Chebyshev{D}}},v::Vector,
             plan=plan_transform(S,v)) = plan*copy(v)


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op){S<:Chebyshev,V<:Chebyshev}(f::ProductFun{S,V}) =
        ProductFun(chebyshevtransform($op(values(f))),space(f))
end



reverseorientation{C<:Chebyshev}(f::Fun{C}) =
    Fun(Chebyshev(reverse(domain(f))),alternatesign!(copy(f.coefficients)))


include("ChebyshevOperators.jl")
