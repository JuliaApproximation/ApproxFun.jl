

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
    @eval ($op){T<:Number,O}(::Type{T},S::UltrasphericalSpace{O})=IFun(($op)(T,1),S)
end


## Transform


#This can be overriden, but the default is to use Chebyshev
transform(sp::DomainSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),ChebyshevSpace(domain(sp)),sp)
itransform(sp::DomainSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,ChebyshevSpace(domain(sp))))


## Algebra

function .*{T,N,a,b}(f::IFun{T,UltrasphericalSpace{a}},g::IFun{N,UltrasphericalSpace{b}})
    @assert domainscompatible(f,g)
    #TODO Coefficient space version
    n = length(f) + length(g) - 1
    f2 = pad(f,n); g2 = pad(g,n)
    
    chop!(IFun(chebyshevtransform(values(f2).*values(g2)),domain(f)),10eps())
end



## Calculus

integrate{T}(f::IFun{T,UltrasphericalSpace{1}})=IFun(fromcanonicalD(f,0)*ultraint(f.coefficients),ChebyshevSpace(domain(f)))




# DirichletSpaces


immutable ChebyshevDirichletSpace{left,right} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end

=={l,r}(a::ChebyshevDirichletSpace{l,r},b::ChebyshevDirichletSpace{l,r})= a.domain==b.domain

function maxspace(a::UltrasphericalSpace,b::ChebyshevDirichletSpace)
    @assert domainscompatible(a,b)
    
    a
end
Base.max(b::ChebyshevDirichletSpace,a::UltrasphericalSpace)=maxspace(a,b)

function minspace(a::UltrasphericalSpace,b::ChebyshevDirichletSpace)
    @assert domainscompatible(a,b)
    
    b
end
minspace(b::ChebyshevDirichletSpace,a::UltrasphericalSpace)=minspace(a,b)