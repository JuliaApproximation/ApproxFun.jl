
export Ultraspherical

#Ultraspherical Spaces



immutable Ultraspherical{O,D<:Domain} <: PolynomialSpace{D}
    domain::D
    Ultraspherical(d)=new(d)
    Ultraspherical()=new(Interval())
end


Base.call{m}(::Type{Ultraspherical{m}})=Ultraspherical{m,Interval{Float64}}()
Base.call{m}(::Type{Ultraspherical{m}},d::Domain)=Ultraspherical{m,typeof(d)}(d)
Base.call{m}(::Type{Ultraspherical{m}},d::Vector)=Ultraspherical{m}(Domain(d))


setdomain{s}(S::Ultraspherical{s},d::Domain)=Ultraspherical{s}(d)


include("Chebyshev.jl")


order{o}(::Ultraspherical{o})=o



## Construction

#domain(S) may be any domain

Base.ones{T<:Number,O}(::Type{T},S::Ultraspherical{O})=Fun(ones(T,1),S)
Base.ones{O}(S::Ultraspherical{O})=Fun(ones(1),S)



## Fast evaluation

Base.first{O,D}(f::Fun{Ultraspherical{O,D}})=foldr(-,coefficients(f,Chebyshev))
Base.last{O,D}(f::Fun{Ultraspherical{O,D}})=reduce(+,coefficients(f,Chebyshev))
identity_fun{m}(d::Ultraspherical{m})=Fun(identity_fun(domain(d)),d)


## Calculus




spacescompatible{O}(a::Ultraspherical{O},b::Ultraspherical{O})=domainscompatible(a,b)
hasfasttransform(::Ultraspherical)=true



include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
include("ContinuousSpace.jl")
