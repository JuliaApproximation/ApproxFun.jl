
export Ultraspherical

#Ultraspherical Spaces



immutable Ultraspherical{O} <: PolynomialSpace
    domain::Union(Interval,AnyDomain)
    Ultraspherical(d)=new(d)
    Ultraspherical()=new(Interval())
end



#Ultraspherical(o::Integer)=Ultraspherical(o,AnyDomain())
#Chebyshev(d::IntervalDomain)=Ultraspherical(0,d)

include("Chebyshev.jl")


order{o}(::Ultraspherical{o})=o



## Construction

#domain(S) may be any domain

Base.ones{T<:Number,O}(::Type{T},S::Ultraspherical{O})=Fun(ones(T,1),S)
Base.ones{O}(S::Ultraspherical{O})=Fun(ones(1),S)



## Fast evaluation

Base.first{O}(f::Fun{Ultraspherical{O}})=foldr(-,coefficients(f,Chebyshev))
Base.last{O}(f::Fun{Ultraspherical{O}})=reduce(+,coefficients(f,Chebyshev))
identity_fun{m}(d::Ultraspherical{m})=Fun(identity_fun(domain(d)),d)


## Calculus

integrate(f::Fun{Ultraspherical{1}})=Fun(fromcanonicalD(f,0)*ultraint(f.coefficients),Chebyshev(domain(f)))



spacescompatible{O}(a::Ultraspherical{O},b::Ultraspherical{O})=domainscompatible(a,b)




include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
include("ContinuousSpace.jl")
