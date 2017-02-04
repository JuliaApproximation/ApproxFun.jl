

export PeriodicInterval

doc"""
	PeriodicInterval(a,b)

represents a periodic interval from `a` to `b`, that is, the point
`b` is identified with `a`.
"""
immutable PeriodicInterval{T} <: PeriodicDomain{T}
    a::T
    b::T
    PeriodicInterval() = new(0,2convert(T,π))
    PeriodicInterval(a,b) = new(a,b)
end

PeriodicInterval()=PeriodicInterval{Float64}()
PeriodicInterval(a::Int,b::Int) = PeriodicInterval(Float64(a),Float64(b)) #convenience method
PeriodicInterval(a,b) = PeriodicInterval{promote_type(typeof(a),typeof(b))}(a,b)
PeriodicInterval(a::Tuple,b::Tuple) = Interval(Vec(a...),Vec(b...))

function Base.convert(::Type{PeriodicInterval},d::ClosedInterval)
	a,b = d.left,d.right
    @assert isfinite(a) && isfinite(b)
    PeriodicInterval(a,b)
end

Segment(d::PeriodicInterval) = Segment(d.a,d.b)
Interval(d::PeriodicInterval) = Interval(d.a,d.b)
PeriodicInterval(d::Segment) = PeriodicInterval(d.a,d.b)

Base.convert{T<:Number}(::Type{PeriodicInterval{T}}, d::PeriodicInterval) = PeriodicInterval{T}(d.a,d.b)

isambiguous(d::PeriodicInterval) = all(isnan(d.a)) && all(isnan(d.b))
Base.convert{T<:Number}(::Type{PeriodicInterval{T}},::AnyDomain) = PeriodicInterval{T}(NaN,NaN)
Base.convert{d,T}(::Type{PeriodicInterval{Vec{d,T}}},::AnyDomain) = PeriodicInterval(Vec(fill(NaN,d)...),Vec(fill(NaN,d)...))
Base.convert{T}(::Type{PeriodicInterval{T}},::AnyDomain) = PeriodicInterval(nan(T),nan(T))


## Information

Base.first(d::PeriodicInterval) = d.a

Base.issubset(a::PeriodicInterval,b::PeriodicInterval) = first(a)∈b && (a.b∈b || a.b == b.b)

# we disable last since the domain is "periodic"


## Map periodic interval


tocanonical{T}(d::PeriodicInterval{T},x) = π*(tocanonical(Segment(d),x)+1)
tocanonicalD{T}(d::PeriodicInterval{T},x) = π*tocanonicalD(Segment(d),x)
fromcanonical(d::PeriodicInterval,θ) = fromcanonical(Segment(d),θ/π-1)
fromcanonicalD(d::PeriodicInterval,θ) = fromcanonicalD(Segment(d),θ/π-1)/π



arclength(d::PeriodicInterval) = norm(d.b - d.a)
Base.angle(d::PeriodicInterval) = angle(d.b - d.a)
complexlength(d::PeriodicInterval) = d.b-d.a
Base.reverse(d::PeriodicInterval) = PeriodicInterval(d.b,d.a)



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b





## algebra

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::PeriodicInterval) = PeriodicInterval($op(c,d.a),$op(c,d.b))
        $op(d::PeriodicInterval,c::Number) = PeriodicInterval($op(d.a,c),$op(d.b,c))
    end
end


@eval /(d::PeriodicInterval,c::Number) = PeriodicInterval(/(d.a,c),/(d.b,c))

+(d1::PeriodicInterval,d2::PeriodicInterval) = PeriodicInterval(d1.a+d2.a,d1.b+d2.b)
