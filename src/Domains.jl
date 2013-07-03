

export Domain,to_uinterval,from_uinterval


abstract Domain


## Interval

type Interval <: Domain
	a
	b
end



function to_uinterval(d::Interval,x)
	(d.a + d.b - 2x)/(d.a - d.b)
end

function to_uintervalD(d::Interval,x)
	 2/( d.b- d.a)
end

function from_uinterval(d::Interval,x)
	.5*(d.a + d.b) + .5*(d.b - d.a)x
end

function from_uintervalD(d::Interval,x)
	 .5*( d.b- d.a)
end

Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b


## Arc

## Polynomial map