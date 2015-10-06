

export PeriodicInterval


immutable PeriodicInterval{T<:Number} <: PeriodicDomain{T}
    a::T
    b::T
    PeriodicInterval()=new(-convert(T,π),convert(T,π))
    PeriodicInterval(a,b)=new(a,b)
end

PeriodicInterval()=PeriodicInterval{Float64}()
PeriodicInterval(a::Int,b::Int) = PeriodicInterval(Float64(a),Float64(b)) #convenience method
PeriodicInterval(a::Number,b::Number) = PeriodicInterval{promote_type(typeof(a),typeof(b))}(a,b)

function PeriodicInterval{T<:Number}(d::AbstractVector{T})
    @assert length(d)==2
    @assert isfinite(d[1]) && isfinite(d[2])
    PeriodicInterval(d...)
end

Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)

Base.convert{T<:Number}(::Type{PeriodicInterval{T}}, d::PeriodicInterval) = PeriodicInterval{T}(d.a,d.b)

isambiguous(d::PeriodicInterval)=isnan(d.a) && isnan(d.b)
Base.convert{T<:Number}(::Type{PeriodicInterval{T}},::AnyDomain)=PeriodicInterval{T}(NaN,NaN)
Base.convert{IT<:PeriodicInterval}(::Type{IT},::AnyDomain)=PeriodicInterval(NaN,NaN)


## Information

Base.first(d::PeriodicInterval)=d.a

Base.issubset(a::PeriodicInterval,b::PeriodicInterval)=first(a)∈b && last(a)∈b

# we disable last since the domain is "periodic"
#Base.last(d::Interval)=d.b


## Map periodic interval


tocanonical{T}(d::PeriodicInterval{T},x)=convert(T,π).*tocanonical(Interval(d),x)
tocanonicalD{T}(d::PeriodicInterval{T},x)=convert(T,π).*tocanonicalD(Interval(d),x)
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
