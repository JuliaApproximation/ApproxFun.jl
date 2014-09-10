

typealias ChebyshevSpace UltrasphericalSpace{0}

## Space conversion default is through Chebyshev

spaceconversion(f::Vector,sp::FunctionSpace)=spaceconversion(f,ChebyshevSpace(AnyDomain()),sp)
spaceconversion(f::Vector,::ChebyshevSpace,sp2::ChebyshevSpace)=f
spaceconversion(f::Vector,::ChebyshevSpace,sp2::FunctionSpace)=error("Override space conversion from ChebyshevSpace to " * string(typeof(sp2)))
spaceconversion(f::Vector,sp2::FunctionSpace,::ChebyshevSpace)=error("Override space conversion from " * string(typeof(sp2)) * " to ChebyshevSpace")
spaceconversion(f::Vector,sp1::FunctionSpace,sp2::FunctionSpace,sp3::FunctionSpace)=spaceconversion(spaceconversion(f,sp1,sp2),sp2,sp3)
spaceconversion(f::Vector,sp1::FunctionSpace,sp2::FunctionSpace)=spaceconversion(f,sp1,ChebyshevSpace(AnyDomain()),sp2)

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

evaluate{T}(f::IFun{T,ChebyshevSpace},x)=clenshaw(f.coefficients,tocanonical(f,x))

## Calculus


# diff T -> U, then convert U -> T
integrate{T}(f::IFun{T,ChebyshevSpace})=IFun(chebyshevintegrate(domain(f),f.coefficients),f.space)
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   


differentiate{T}(f::IFun{T,ChebyshevSpace})=IFun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(IFun(x->tocanonicalD(d,x),d).*IFun(diff(IFun(cfs)),d)).coefficients