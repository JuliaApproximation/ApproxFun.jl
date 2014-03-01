

export Interval



## Standard interval

type Interval{T<:Number} <: IntervalDomain
	a::T
	b::T
end

Interval()=Interval(-1.,1.)


Base.convert{D<:Interval}(::Type{D},i::Vector)=Interval(i)
Interval(a::Number,b::Number) = Interval(promote(a,b)...)


## Information

Base.first(d::Interval)=d.a
Base.last(d::Interval)=d.b



## Map interval



tocanonical(d::Interval,x)=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD(d::Interval,x)=2/( d.b- d.a)
fromcanonical(d::Interval,x)=.5*(d.a + d.b) + .5*(d.b - d.a)x
fromcanonicalD(d::Interval,x)=.5*( d.b- d.a)



Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b



##Integration and differentiation


# diff T -> U, then convert U -> T
Base.diff{T<:Number,M<:Interval}(f::IFun{T,M})=tocanonicalD(f.domain,0)*IFun(ultraiconversion(ultradiff(f.coefficients)),f.domain)
integrate{T<:Number,M<:Interval}(f::IFun{T,M})=fromcanonicalD(f.domain,0)*IFun(ultraint(ultraconversion(f.coefficients)),f.domain)    


Base.diff{T<:Number,M<:IntervalDomain}(f::IFun{T,M})=IFun(x->tocanonicalD(f.domain,x),f.domain).*IFun(diff(IFun(f)),f.domain)

function Base.diff(f::IFun,k::Integer)
    @assert k >= 0
    (k==0)?f:diff(diff(f),k-1)
end