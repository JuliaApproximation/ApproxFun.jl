

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

##max space



function maxspace{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder > border?a:b
end

function minspace{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end



## Construction

#domain(S) may be any domain
for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number,O}(::Type{T},S::UltrasphericalSpace{O})=IFun(($op)(T,1),S)
end


## Transform


#This can be overriden, but the default is to use Chebyshev
transform(sp::DomainSpace,vals::Vector)=spaceconversion(chebyshevtransform(vals),ChebyshevSpace(domain(sp)),sp)
itransform(sp::DomainSpace,cfs::Vector)=ichebyshevtransform(spaceconversion(cfs,sp,ChebyshevSpace(domain(sp))))




## Calculus

integrate{T}(f::IFun{T,UltrasphericalSpace{1}})=IFun(fromcanonicalD(f,0)*ultraint(f.coefficients),ChebyshevSpace(domain(f)))




# DirichletSpaces


immutable ChebyshevDirichletSpace{T<:Union(IntervalDomain,AnyDomain)} <: IntervalDomainSpace
    domain::T 
    left::Int
    right::Int    
end

==(a::ChebyshevDirichletSpace,b::ChebyshevDirichletSpace)= a.domain==b.domain && a.left==b.left && a.right==b.right

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