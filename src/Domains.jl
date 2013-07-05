

export Domain,tointerval,frominterval


abstract Domain


## Interval

type Interval{T<:Real} <: Domain
	a::T
	b::T
end

tointerval(d::Interval,x)=(d.a + d.b - 2x)/(d.a - d.b)
tointervalD(d::Interval,x)=2/( d.b- d.a)
frominterval(d::Interval,x)=.5*(d.a + d.b) + .5*(d.b - d.a)x
fromintervalD(d::Interval,x)=.5*( d.b- d.a)

for op = (:tointerval,:tointervalD,:frominterval,:fromintervalD)
    @eval ($op)(f::AbstractFun,x)=($op)(f.domain,x)
end

Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b


## Arc

## Polynomial map