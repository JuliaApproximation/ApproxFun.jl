
export UltrasphericalSpace

#Ultraspherical Spaces



immutable UltrasphericalSpace{O} <: IntervalSpace
    domain::Union(Interval,AnyDomain)
    UltrasphericalSpace(d)=new(d)
    UltrasphericalSpace()=new(Interval())
end



#UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,AnyDomain())
#Chebyshev(d::IntervalDomain)=UltrasphericalSpace(0,d)

include("Chebyshev.jl")


order{o}(::UltrasphericalSpace{o})=o



## Construction

#domain(S) may be any domain

Base.ones{T<:Number,O}(::Type{T},S::UltrasphericalSpace{O})=Fun(ones(T,1),S)
Base.ones{O}(S::UltrasphericalSpace{O})=Fun(ones(1),S)    


## Transform


#This can be overriden, but the default is to use Chebyshev
# transform(sp::IntervalSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),Chebyshev(domain(sp)),sp)
# itransform(sp::IntervalSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,Chebyshev(domain(sp))))



## Fast evaluation

Base.first{O}(f::Fun{UltrasphericalSpace{O}})=foldr(-,canonicalcoefficients(f))
Base.last{O}(f::Fun{UltrasphericalSpace{O}})=reduce(+,canonicalcoefficients(f))
identity_fun{m}(d::UltrasphericalSpace{m})=Fun(identity_fun(domain(d)),d)


## Calculus

integrate(f::Fun{UltrasphericalSpace{1}})=Fun(fromcanonicalD(f,0)*ultraint(f.coefficients),Chebyshev(domain(f)))







include("UltrasphericalOperators.jl")
include("DirichletSpace.jl")
