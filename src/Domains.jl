

export Domain,to_uinterval,from_uinterval


abstract Domain


type Interval <: Domain
	a
	b
end


function to_uinterval(d::Interval,x)
	(d.a + d.b - 2x)/(d.a - d.b)
end

function from_uinterval(d::Interval,x)
	.5*(d.a + d.b) + .5*(d.b - d.a)x
end