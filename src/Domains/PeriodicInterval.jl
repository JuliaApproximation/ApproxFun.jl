

export PeriodicInterval


immutable PeriodicInterval{T<:Number} <: PeriodicDomain{T}
    a::T
    b::T
end

PeriodicInterval()=PeriodicInterval(-1.π,1.π)
PeriodicInterval(a::Int,b::Int) = PeriodicInterval(float64(a),float64(b))   #convenience method

Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)

Base.convert{T<:Number}(::Type{PeriodicInterval{T}}, d::PeriodicInterval) = PeriodicInterval{T}(d.a,d.b)

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



Base.length(d::PeriodicInterval) = abs(d.b - d.a)
Base.angle(d::PeriodicInterval) = angle(d.b - d.a)
Base.reverse(d::PeriodicInterval)=PeriodicInterval(d.b,d.a)



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b





## algebra

for op in (:*,:+,:-,:.*,:.+,:.-)
    @eval begin
        $op(c::Number,d::PeriodicInterval)=PeriodicInterval($op(c,d.a),$op(c,d.b))
        $op(d::PeriodicInterval,c::Number)=PeriodicInterval($op(d.a,c),$op(d.b,c))
    end
end

for op in (:/,:./)
    @eval $op(d::PeriodicInterval,c::Number)=PeriodicInterval($op(d.a,c),$op(d.b,c))
end


+(d1::PeriodicInterval,d2::PeriodicInterval)=PeriodicInterval(d1.a+d2.a,d1.b+d2.b)

