
export Ultraspherical

#Ultraspherical Spaces



immutable Ultraspherical{O,D<:Domain} <: PolynomialSpace{D}
    domain::D
    Ultraspherical(d) = new(d)
    Ultraspherical() = new(Interval())
end


@compat (::Type{Ultraspherical{m}}){m}() = (@assert m ≠ 0;Ultraspherical{m,Interval{Float64}}())
@compat (::Type{Ultraspherical{m}}){m}(d::Domain) = (@assert m ≠ 0;Ultraspherical{m,typeof(d)}(d))
@compat (::Type{Ultraspherical{m}}){m}(d::Vector) = (@assert m ≠ 0;Ultraspherical{m}(Domain(d)))


setdomain{s}(S::Ultraspherical{s},d::Domain) = Ultraspherical{s}(d)


order{o}(::Ultraspherical{o}) = o
order{o,D}(::Type{Ultraspherical{o,D}}) = o


canonicalspace(S::Ultraspherical) = Chebyshev(domain(S))


## Construction

#domain(S) may be any domain

Base.ones{T<:Number,O}(::Type{T},S::Ultraspherical{O}) = Fun(ones(T,1),S)
Base.ones{O}(S::Ultraspherical{O}) = Fun(ones(1),S)



## Fast evaluation

Base.first{O,D}(f::Fun{Ultraspherical{O,D}}) = isinteger(O)         ?
                                foldr(-,coefficients(f,Chebyshev))  :
                                f(first(domain(f)))
Base.last{O,D}(f::Fun{Ultraspherical{O,D}}) = isinteger(O) ?
                                reduce(+,coefficients(f,Chebyshev)) :
                                f(last(domain(f)))
identity_fun{m}(d::Ultraspherical{m}) = Fun(identity_fun(domain(d)),d)


## Calculus




spacescompatible{O}(a::Ultraspherical{O},b::Ultraspherical{O}) = domainscompatible(a,b)
hasfasttransform(::Ultraspherical) = true



include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
include("ContinuousSpace.jl")
