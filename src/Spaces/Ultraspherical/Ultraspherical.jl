
export Ultraspherical

#Ultraspherical Spaces



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


## Construction

#domain(S) may be any domain

Base.ones{T<:Number}(::Type{T},S::Ultraspherical) = Fun(ones(T,1),S)
Base.ones(S::Ultraspherical) = Fun(ones(1),S)



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
