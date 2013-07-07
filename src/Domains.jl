

export Domain,tocanonical,fromcanonical,Circle


abstract Domain


##General routines

for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::AbstractFun,x)=($op)(f.domain,x)
end






## Interval Domains

abstract IntervalDomain <: Domain

chebyshevpoints(n::Integer)= cos(1.π*[n-1:-1:0]/(n-1))
points(d::IntervalDomain,n::Integer) = fromcanonical(d, chebyshevpoints(n))


type Interval{T<:Number} <: IntervalDomain
	a::T
	b::T
end

Interval()=Interval(-1.,1.)


## Points








## Map interval



tocanonical(d::Interval,x)=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD(d::Interval,x)=2/( d.b- d.a)
fromcanonical(d::Interval,x)=.5*(d.a + d.b) + .5*(d.b - d.a)x
fromcanonicalD(d::Interval,x)=.5*( d.b- d.a)



Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b



## Arc

## Polynomial map



####
## Periodic domains

abstract PeriodicDomain <: Domain

points(d::PeriodicDomain,n::Integer) = fromcanonical(d, fourierpoints(n))


type PeriodicInterval{T<:Number} <: PeriodicDomain
	a::T
	b::T
end

PeriodicInterval()=PeriodicInterval(-1.π,1.π)


Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)

fourierpoints(n::Integer)= 1.π*[-1.:2/n:1. - 2/n]





## Map periodic interval


tocanonical(d::PeriodicInterval,x)=1.π.*tocanonical(Interval(d),x)
tocanonicalD(d::PeriodicInterval,x)=1.π.*tocanonicalD(Interval(d),x)
fromcanonical(d::PeriodicInterval,θ)=fromcanonical(Interval(d),θ/π)
fromcanonicalD(d::PeriodicInterval,θ)=fromcanonicalD(Interval(d),θ/π)/π



Base.length(d::PeriodicInterval) = d.b - d.a



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b


##  Circle


type Circle <: PeriodicDomain
	center
	radius
end

Circle(r)=Circle(0.,r)
Circle()=Circle(1.)


tocanonical(d::Circle,ζ)=-1.im*log((ζ-d.center)/d.radius)
tocanonicalD(d::Circle,ζ)=-1.im./((ζ-d.center))  #TODO: Check formula
fromcanonical(d::Circle,θ)=d.radius*exp(1.im*θ) + d.center
fromcanonicalD(d::Circle,θ)=d.radius*1.im*exp(1.im*θ)



Base.length(d::Circle) = 2π*d.radius



==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius

