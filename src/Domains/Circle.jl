

export Circle


##  Circle


immutable Circle{T<:Number,V<:Real} <: PeriodicDomain{Complex{Float64}}
	center::T
	radius::V
end

Circle(r)=Circle(0.,r)
Circle()=Circle(1.)

function Circle{T<:Number,V<:Real}(c::Vector{T},r::Vector{V})
    @assert length(c) == length(r)
    if length(c) == 1
        Circle(c[1],r[1])
    else
        [Circle(c[1],r[1]);Circle(c[2:end],r[2:end])]
    end
end
Circle{T<:Number,V<:Real}(c::Vector{T},r::V) = Circle(c,ones(V,length(c))r)

function tocanonical(d::Circle,ζ)
    v=mappoint(d,Circle(),ζ)- 0.im#Subtract 0.im so branch cut is right

    -1.im.*log(v)
end

tocanonicalD(d::Circle,ζ)=-1.im./(ζ.-d.center)  #TODO: Check formula
fromcanonical(d::Circle,θ)=d.radius*exp(1.im*θ) .+ d.center
fromcanonicalD(d::Circle,θ)=d.radius*1.im*exp(1.im*θ)

canonicaldomain(d::Circle)=PeriodicInterval()

Base.in(z,d::Circle)=isapprox(abs(z-d.center),d.radius)

Base.length(d::Circle) = 2π*d.radius



==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius



function mappoint(d1::Circle,d2::Circle,z)
   v=(z-d1.center)/d1.radius
   v*d2.radius+d2.center
end



