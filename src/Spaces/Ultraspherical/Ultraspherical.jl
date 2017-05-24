
export Ultraspherical

#Ultraspherical Spaces


doc"""
`Ultraspherical(λ)` is the space spanned by the ultraspherical polynomials
```
    C_0^{(λ)}(x),C_1^{(λ)}(x),C_2^{(λ)}(x),…
```
Note that `λ=1` this reduces to Chebyshev polynomials of the second kind:
`C_k^{(1)}(x) = U_k(x)`.
"""
struct Ultraspherical{T,D<:Domain,R} <: PolynomialSpace{D,R}
    order::T
    domain::D
    Ultraspherical{T,D,R}(m::T,d::D) where {T,D,R} = (@assert m ≠ 0; new(m,d))
    Ultraspherical{T,D,R}(m::Number,d::Domain) where {T,D,R} = (@assert m ≠ 0; new(convert(T,m),convert(D,d)))
    Ultraspherical{T,D,R}(d::Domain) where {T,D,R} = new(one(T),convert(D,d))
    Ultraspherical{T,D,R}(m::Number) where {T,D,R} = (@assert m ≠ 0; new(T(m),D()))
end

Ultraspherical(m::Number,d::Domain) = Ultraspherical{typeof(m),typeof(d),real(prectype(d))}(m,d)
Ultraspherical(m::Number,d) = Ultraspherical(m,Domain(d))
Ultraspherical(m::Number) = Ultraspherical(m,Interval())


order(S::Ultraspherical) = S.order
setdomain(S::Ultraspherical,d::Domain) = Ultraspherical(order(S),d)





canonicalspace(S::Ultraspherical) = Chebyshev(domain(S))


struct UltrasphericalPlan{CT,FT}
    chebplan::CT
    cheb2legplan::FT

    (::Type{UltrasphericalPlan{CT,FT}}){CT,FT}(cp,c2lp) = new{CT,FT}(cp,c2lp)
end

struct UltrasphericalIPlan{CT,FT}
    chebiplan::CT
    leg2chebplan::FT

    (::Type{UltrasphericalIPlan{CT,FT}}){CT,FT}(cp,c2lp) = new{CT,FT}(cp,c2lp)
end

function UltrasphericalPlan(λ::Number,vals)
    if λ == 0.5
        cp = plan_transform(Chebyshev(),vals)
        c2lp = FastTransforms.th_cheb2legplan(eltype(vals),length(vals))
        UltrasphericalPlan{typeof(cp),typeof(c2lp)}(cp,c2lp)
    else
        error("Not implemented")
    end
end

function UltrasphericalIPlan(λ::Number,cfs)
    if λ == 0.5
        cp=plan_itransform(Chebyshev(),cfs)
        c2lp=FastTransforms.th_leg2chebplan(eltype(cfs),length(cfs))
        UltrasphericalIPlan{typeof(cp),typeof(c2lp)}(cp,c2lp)
    else
        error("Not implemented")
    end
end

*(UP::UltrasphericalPlan,v::AbstractVector) =
    UP.cheb2legplan*(UP.chebplan*v)
*(UP::UltrasphericalIPlan,v::AbstractVector) =
    UP.chebiplan*(UP.leg2chebplan*v)


plan_transform(sp::Ultraspherical{Int},vals::Vector) = CanonicalTransformPlan(sp,vals)
plan_transform(sp::Ultraspherical,vals::Vector) = UltrasphericalPlan(order(sp),vals)
plan_itransform(sp::Ultraspherical{Int},cfs::Vector) = ICanonicalTransformPlan(sp,cfs)
plan_itransform(sp::Ultraspherical,cfs::Vector) = UltrasphericalIPlan(order(sp),cfs)

## Construction

#domain(S) may be any domain

Base.ones{T<:Number}(::Type{T},S::Ultraspherical) = Fun(S,ones(T,1))
Base.ones(S::Ultraspherical) = Fun(S,ones(1))



## Fast evaluation

function Base.first{D,R}(f::Fun{Ultraspherical{Int,D,R}})
    n = length(f.coefficients)
    n == 0 && return zero(eltype(f))
    n == 1 && return first(f.coefficients)
    foldr(-,coefficients(f,Chebyshev))
end

Base.last{D,R}(f::Fun{Ultraspherical{Int,D,R}}) = reduce(+,coefficients(f,Chebyshev))

Base.first{O,D,R}(f::Fun{Ultraspherical{O,D,R}}) = f(first(domain(f)))
Base.last{O,D,R}(f::Fun{Ultraspherical{O,D,R}}) = f(last(domain(f)))

identity_fun(d::Ultraspherical) = Fun(identity_fun(domain(d)),d)


## Calculus




spacescompatible(a::Ultraspherical,b::Ultraspherical) =
    order(a) == order(b) && domainscompatible(a,b)
hasfasttransform(::Ultraspherical) = true



include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
include("ContinuousSpace.jl")
