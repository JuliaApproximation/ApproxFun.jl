

export PeriodicInterval


immutable PeriodicInterval{T<:Number} <: PeriodicDomain
	a::T
	b::T
end

PeriodicInterval()=PeriodicInterval(-1.π,1.π)


Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)



function PeriodicInterval{T<:Number}(d::Vector{T})
    @assert length(d) == 2
    
    if abs(d[1]) ==Inf
        PeriodicLine(d)
    else
        PeriodicInterval(d[1],d[2])
    end
end

## Information

Base.first(d::PeriodicInterval)=d.a

# we disable last since the domain is "periodic"
#Base.last(d::Interval)=d.b


## Map periodic interval


tocanonical(d::PeriodicInterval,x)=1.π.*tocanonical(Interval(d),x)
tocanonicalD(d::PeriodicInterval,x)=1.π.*tocanonicalD(Interval(d),x)
fromcanonical(d::PeriodicInterval,θ)=fromcanonical(Interval(d),θ/π)
fromcanonicalD(d::PeriodicInterval,θ)=fromcanonicalD(Interval(d),θ/π)/π



Base.length(d::PeriodicInterval) = d.b - d.a



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b



