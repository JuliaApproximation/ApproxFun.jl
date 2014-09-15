

export Circle


##  Circle


immutable Circle{T<:Number} <: PeriodicDomain
	center::T
	radius::Float64
end

Circle(r)=Circle(0.,r)
Circle()=Circle(1.)


function tocanonical(d::Circle,ζ)
    v=(ζ-d.center)/d.radius - 0.im#Subtract 0.im so branch cut is right
    -1.im.*log(v)
end

tocanonicalD(d::Circle,ζ)=-1.im./(ζ.-d.center)  #TODO: Check formula
fromcanonical(d::Circle,θ)=d.radius*exp(1.im*θ) .+ d.center
fromcanonicalD(d::Circle,θ)=d.radius*1.im*exp(1.im*θ)



Base.length(d::Circle) = 2π*d.radius



==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius






