
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


## Transform


#This can be overriden, but the default is to use Chebyshev
# transform(sp::IntervalSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),Chebyshev(domain(sp)),sp)
# itransform(sp::IntervalSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,Chebyshev(domain(sp))))



## Fast evaluation

Base.first{O}(f::Fun{Ultraspherical{O}})=foldr(-,canonicalcoefficients(f))
Base.last{O}(f::Fun{Ultraspherical{O}})=reduce(+,canonicalcoefficients(f))
identity_fun{m}(d::Ultraspherical{m})=Fun(identity_fun(domain(d)),d)


## Calculus

integrate(f::Fun{Ultraspherical{1}})=Fun(fromcanonicalD(f,0)*ultraint(f.coefficients),Chebyshev(domain(f)))







include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
