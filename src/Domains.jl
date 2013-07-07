

export Domain,tocanonical,fromcanonical


abstract Domain





## Interval



type Interval{T<:Number} <: Domain
	a::T
	b::T
end

type PeriodicInterval{T<:Number} <: Domain
	a::T
	b::T
end

Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)


## Points

chebyshevpoints(n::Integer)= cos(1.π*[n-1:-1:0]/(n-1))
fourierpoints(n::Integer)= 1.π*[-1.:2/n:1. - 2/n]

points(d::Interval,n::Integer) = fromcanonical(d, chebyshevpoints(n))
points(d::PeriodicInterval,n::Integer) = fromcanonical(d, fourierpoints(n))


##Fun calls to domain calls

for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::AbstractFun,x)=($op)(f.domain,x)
end


## Map interval



tocanonical(d::Interval,x)=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD(d::Interval,x)=2/( d.b- d.a)
fromcanonical(d::Interval,x)=.5*(d.a + d.b) + .5*(d.b - d.a)x
fromcanonicalD(d::Interval,x)=.5*( d.b- d.a)



Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b


## Map periodic interval


tocanonical(d::PeriodicInterval,x)=1.π.*tocanonical(Interval(d),x)
tocanonicalD(d::PeriodicInterval,x)=1.π.*tocanonicalD(Interval(d),x)
fromcanonical(d::PeriodicInterval,x)=fromcanonical(Interval(d),x/π)
fromcanonicalD(d::PeriodicInterval,x)=fromcanonicalD(Interval(d),x/π)/π



Base.length(d::PeriodicInterval) = d.b - d.a



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b


## Arc

## Polynomial map




