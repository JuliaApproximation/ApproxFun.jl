

#Ultraspherical Spaces



type UltrasphericalSpace{T<:IntervalDomain} <: IntervalDomainSpace
    order::Int
    domain::T     
end

#UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,AnyDomain())
ChebyshevSpace(d::IntervalDomain)=UltrasphericalSpace(0,d)




#TODO: bad override?
==(a::UltrasphericalSpace,b::UltrasphericalSpace)=a.order==b.order && domainscompatible(a,b)

##max space



function maxspace(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order > b.order?a:b
end

function minspace(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order < b.order?a:b
end



## Operator space manipulation



# DirichletSpaces


type ChebyshevDirichletSpace{T<:Union(IntervalDomain,AnyDomain)} <: IntervalDomainSpace
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


##Integration and differentiation


# diff T -> U, then convert U -> T
function differentiate(sp::UltrasphericalSpace,cfs::Vector)
    @assert sp.order==0
    chebyshevdifferentiate(domain(sp),cfs)
end

function integrate(sp::UltrasphericalSpace,cfs::Vector)
    @assert sp.order==0
    chebyshevintegrate(domain(sp),cfs)
end

function chebyshevsum(d::IntervalDomain,cfs::Vector)
    cf=IFun(chebyshevintegrate(d,cfs),d)
    last(cf) - first(cf)
end


function Base.sum(sp::UltrasphericalSpace,cfs::Vector)
    @assert sp.order==0
    chebyshevsum(domain(sp),cfs)
end

chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
##TODO: can we get rid of reference to IFun?
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(IFun(x->tocanonicalD(d,x),d).*IFun(diff(IFun(cfs)),d)).coefficients
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   



