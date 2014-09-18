
export UltrasphericalSpace

#Ultraspherical Spaces



immutable UltrasphericalSpace{O} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end

#UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,AnyDomain())
#ChebyshevSpace(d::IntervalDomain)=UltrasphericalSpace(0,d)

include("ChebyshevSpace.jl")


order{o}(::UltrasphericalSpace{o})=o



## Construction

#domain(S) may be any domain
for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number,O}(::Type{T},S::UltrasphericalSpace{O})=Fun(($op)(T,1),S)
    @eval ($op){O}(S::UltrasphericalSpace{O})=Fun(($op)(1),S)    
end


## Transform


#This can be overriden, but the default is to use Chebyshev
transform(sp::IntervalDomainSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),ChebyshevSpace(domain(sp)),sp)
itransform(sp::IntervalDomainSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,ChebyshevSpace(domain(sp))))



## Fast evaluation

Base.first{O}(f::Fun{UltrasphericalSpace{O}})=foldr(-,canonicalcoefficients(f))
Base.last{O}(f::Fun{UltrasphericalSpace{O}})=reduce(+,canonicalcoefficients(f))



## Calculus

integrate(f::Fun{UltrasphericalSpace{1}})=Fun(fromcanonicalD(f,0)*ultraint(f.coefficients),ChebyshevSpace(domain(f)))




# DirichletSpaces


immutable ChebyshevDirichletSpace{left,right} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end






include("UltrasphericalOperators.jl")