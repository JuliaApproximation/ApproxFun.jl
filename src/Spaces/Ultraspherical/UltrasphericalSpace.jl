
export UltrasphericalSpace

#Ultraspherical Spaces



immutable UltrasphericalSpace{O} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end

#UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,AnyDomain())
#ChebyshevSpace(d::IntervalDomain)=UltrasphericalSpace(0,d)

include("ChebyshevSpace.jl")


order{o}(::UltrasphericalSpace{o})=o



#TODO: bad override?
=={T}(a::UltrasphericalSpace{T},b::UltrasphericalSpace{T})=domainscompatible(a,b)




## Construction

#domain(S) may be any domain
for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number,O}(::Type{T},S::UltrasphericalSpace{O})=Fun(($op)(T,1),S)
end


## Transform


#This can be overriden, but the default is to use Chebyshev
transform(sp::IntervalDomainSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),ChebyshevSpace(domain(sp)),sp)
itransform(sp::IntervalDomainSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,ChebyshevSpace(domain(sp))))



## Fast evaluation

Base.first{T,O}(f::Fun{T,UltrasphericalSpace{O}})=foldr(-,coefficients(f))
Base.last{T,O}(f::Fun{T,UltrasphericalSpace{O}})=reduce(+,coefficients(f))



## Calculus

integrate{T}(f::Fun{T,UltrasphericalSpace{1}})=Fun(fromcanonicalD(f,0)*ultraint(f.coefficients),ChebyshevSpace(domain(f)))




# DirichletSpaces


immutable ChebyshevDirichletSpace{left,right} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end

=={l,r}(a::ChebyshevDirichletSpace{l,r},b::ChebyshevDirichletSpace{l,r})= a.domain==b.domain




include("UltrasphericalOperators.jl")