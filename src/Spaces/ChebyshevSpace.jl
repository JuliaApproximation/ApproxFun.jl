

typealias ChebyshevSpace UltrasphericalSpace{0}



## Space conversion default is through Chebyshev

spaceconversion(f::Vector,sp::FunctionSpace)=spaceconversion(f,ChebyshevSpace(AnyDomain()),sp)
spaceconversion{T<:FunctionSpace}(f::Vector,::T,sp2::T)=f       #TODO: Check domains?
spaceconversion(f::Vector,sp1::FunctionSpace,sp2::FunctionSpace,sp3::FunctionSpace)=spaceconversion(spaceconversion(f,sp1,sp2),sp2,sp3)


## spaceconversion defaults to calling ConversionOperator, otherwise it tries to pipe through ChebyshevSpace

function spaceconversion{A<:FunctionSpace,B<:FunctionSpace}(f::Vector,a::A,b::B)
    ct=conversion_type(a,b)
    if ct==a
        ConversionOperator(a,b)*f
    elseif ct==b
        ConversionOperator(b,a)\f    
    elseif typeof(a) <: ChebyshevSpace
        error("Override spaceconversion or implement ConversionOperator from ChebyshevSpace to " * string(B))
    elseif typeof(b) <: ChebyshevSpace
        error("Override spaceconversion or implement ConversionOperator from " * string(A) * " to ChebyshevSpace")
    else
        spaceconversion(f,a,ChebyshevSpace(AnyDomain()),b)
    end
end

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