
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
immutable Ultraspherical{T,D<:Domain} <: PolynomialSpace{D}
    order::T
    domain::D
    Ultraspherical(m::T,d::D) = (@assert m ≠ 0; new(m,d))
    Ultraspherical(m::Number,d::Domain) = (@assert m ≠ 0; new(convert(T,m),convert(D,d)))
    Ultraspherical(d::Domain) = new(one(T),convert(D,d))
    Ultraspherical(m::Number) = (@assert m ≠ 0; new(T(m),D()))
end

Ultraspherical(m::Number,d::Domain) = Ultraspherical{typeof(m),typeof(d)}(m,d)
Ultraspherical(m::Number,d) = Ultraspherical(m,Domain(d))
Ultraspherical(m::Number) = Ultraspherical(m,Interval())


order(S::Ultraspherical) = S.order
setdomain(S::Ultraspherical,d::Domain) = Ultraspherical(order(S),d)





canonicalspace(S::Ultraspherical) = Chebyshev(domain(S))


immutable UltrasphericalPlan{CT,FT}
    chebplan::CT
    cheb2legplan::FT

    UltrasphericalPlan(cp,c2lp) = new(cp,c2lp)
end

immutable UltrasphericalIPlan{CT,FT}
    chebiplan::CT
    leg2chebplan::FT

    UltrasphericalIPlan(cp,c2lp) = new(cp,c2lp)
end

function UltrasphericalPlan(λ::Number,vals)
    if λ == 0.5
        cp=plan_transform(Chebyshev(),vals)
        c2lp=FastTransforms.th_cheb2legplan(eltype(vals),length(vals))
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
    UP.cheb2legplan*transform(Chebyshev(),v,UP.chebplan)
*(UP::UltrasphericalIPlan,v::AbstractVector) =
    itransform(Chebyshev(),UP.leg2chebplan*v,UP.chebiplan)


plan_transform(S::Ultraspherical{Int},vals::Vector) = plan_transform(canonicalspace(S),vals)
plan_transform(S::Ultraspherical,vals::Vector) = UltrasphericalPlan(order(S),vals)
plan_itransform(S::Ultraspherical{Int},cfs::Vector) = plan_itransform(canonicalspace(S),cfs)
plan_itransform(S::Ultraspherical,cfs::Vector) = UltrasphericalIPlan(order(S),cfs)

transform(S::Ultraspherical,vals,pl::UltrasphericalPlan) = pl*vals
itransform(S::Ultraspherical,cfs,pl::UltrasphericalIPlan) = pl*cfs

## Construction

#domain(S) may be any domain

Base.ones{T<:Number}(::Type{T},S::Ultraspherical) = Fun(S,ones(T,1))
Base.ones(S::Ultraspherical) = Fun(S,ones(1))



## Fast evaluation

Base.first{D}(f::Fun{Ultraspherical{Int,D}}) = foldr(-,coefficients(f,Chebyshev))
Base.last{D}(f::Fun{Ultraspherical{Int,D}}) = reduce(+,coefficients(f,Chebyshev))

Base.first{O,D}(f::Fun{Ultraspherical{O,D}}) = f(first(domain(f)))
Base.last{O,D}(f::Fun{Ultraspherical{O,D}}) = f(last(domain(f)))

identity_fun(d::Ultraspherical) = Fun(identity_fun(domain(d)),d)


## Calculus




spacescompatible(a::Ultraspherical,b::Ultraspherical) =
    order(a) == order(b) && domainscompatible(a,b)
hasfasttransform(::Ultraspherical) = true



include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
include("ContinuousSpace.jl")
