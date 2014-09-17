

typealias ChebyshevSpace UltrasphericalSpace{0}


Space(d::IntervalDomain)=ChebyshevSpace(d)
canonicalspace(S::IntervalDomainSpace)=ChebyshevSpace(domain(S))

function spaceconversion(g::Vector,::ConstantSpace,::ChebyshevSpace)
    @assert length(g)==1
    g
end

function spaceconversion(g::Vector,::ChebyshevSpace,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::ChebyshevSpace,vals::Vector)=chebyshevtransform(vals)
itransform(::ChebyshevSpace,cfs::Vector)=ichebyshevtransform(cfs)


## Evaluation

evaluate{T}(f::Fun{T,ChebyshevSpace},x)=clenshaw(f.coefficients,tocanonical(f,x))

## Calculus


# diff T -> U, then convert U -> T
integrate{T}(f::Fun{T,ChebyshevSpace})=Fun(chebyshevintegrate(domain(f),f.coefficients),f.space)
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   


differentiate{T}(f::Fun{T,ChebyshevSpace})=Fun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(Fun(x->tocanonicalD(d,x),d).*Fun(diff(Fun(cfs)),d)).coefficients


## identity_fun

identity_fun(d::ChebyshevSpace)=identity_fun(domain(d))