

#Ultraspherical Spaces



immutable UltrasphericalSpace{O} <: IntervalDomainSpace
    domain::IntervalDomain
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



## Operator space manipulation



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


## evaluation

evaluate(sp::ChebyshevSpace,cfs::Vector,x)=clenshaw(cfs,tocanonical(sp,x))


##Integration and differentiation


# diff T -> U, then convert U -> T
function differentiate(sp::ChebyshevSpace,cfs::Vector)
    chebyshevdifferentiate(domain(sp),cfs)
end

integrate(sp::ChebyshevSpace,cfs::Vector)=chebyshevintegrate(domain(sp),cfs)
integrate(sp::UltrasphericalSpace{1},cfs::Vector)=fromcanonicalD(domain(sp),0)*ultraint(cfs)

function chebyshevsum(d::IntervalDomain,cfs::Vector)
    cf=IFun(chebyshevintegrate(d,cfs),d)
    last(cf) - first(cf)
end


function Base.sum(sp::ChebyshevSpace,cfs::Vector)
    chebyshevsum(domain(sp),cfs)
end

chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
##TODO: can we get rid of reference to IFun?
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(IFun(x->tocanonicalD(d,x),d).*IFun(diff(IFun(cfs)),d)).coefficients
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   





