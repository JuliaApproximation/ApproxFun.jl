

#Ultraspherical Spaces



immutable UltrasphericalSpace{O} <: IntervalDomainSpace
    domain::Union(IntervalDomain,AnyDomain)
end

#UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,AnyDomain())
#ChebyshevSpace(d::IntervalDomain)=UltrasphericalSpace(0,d)


typealias ChebyshevSpace UltrasphericalSpace{0}


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


## transform

itransform(::ChebyshevSpace,cfs::Vector)=ichebyshevtransform(cfs)
itransform(sp::UltrasphericalSpace,cfs::Vector)=itransform(ChebyshevSpace(domain(sp)),spaceconversion(cfs,sp,ChebyshevSpace(domain(sp))))


## evaluation

evaluate{T}(f::IFun{T,ChebyshevSpace},x)=clenshaw(f.coefficients,tocanonical(f,x))


##Integration and differentiation


# diff T -> U, then convert U -> T
integrate{T}(f::IFun{T,ChebyshevSpace})=IFun(chebyshevintegrate(domain(f),f.coefficients),f.space)
integrate{T}(f::IFun{T,UltrasphericalSpace{1}})=IFun(fromcanonicalD(f,0)*ultraint(f.coefficients),ChebyshevSpace(domain(f)))
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   


differentiate{T}(f::IFun{T,ChebyshevSpace})=IFun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(IFun(x->tocanonicalD(d,x),d).*IFun(diff(IFun(cfs)),d)).coefficients










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