

export PeriodicInterval

doc"""
	PeriodicInterval(a,b)

represents a periodic interval from `a` to `b`, that is, the point
`b` is identified with `a`.
"""
struct PeriodicInterval{T} <: PeriodicDomain{T}
    a::T
    b::T
    PeriodicInterval{T}() where {T} = new{T}(0,2convert(T,π))
    PeriodicInterval{T}(a,b) where {T} = new{T}(a,b)
end

PeriodicInterval()=PeriodicInterval{Float64}()
PeriodicInterval(a::Int,b::Int) = PeriodicInterval(Float64(a),Float64(b)) #convenience method
PeriodicInterval(a,b) = PeriodicInterval{promote_type(typeof(a),typeof(b))}(a,b)
PeriodicInterval(a::Tuple,b::Tuple) = Interval(Vec(a...),Vec(b...))

function convert(::Type{PeriodicInterval},d::ClosedInterval)
	a,b = d.left,d.right
    @assert isfinite(a) && isfinite(b)
    PeriodicInterval(a,b)
end

Segment(d::PeriodicInterval) = Segment(d.a,d.b)
Interval(d::PeriodicInterval) = Interval(d.a,d.b)
PeriodicInterval(d::Segment) = PeriodicInterval(d.a,d.b)

convert(::Type{PeriodicInterval{T}}, d::PeriodicInterval) where {T<:Number} = PeriodicInterval{T}(d.a,d.b)

isambiguous(d::PeriodicInterval) = all(isnan(d.a)) && all(isnan(d.b))
convert(::Type{PeriodicInterval{T}},::AnyDomain) where {T<:Number} = PeriodicInterval{T}(NaN,NaN)
convert(::Type{PeriodicInterval{Vec{d,T}}},::AnyDomain) where {d,T} = PeriodicInterval(Vec(fill(NaN,d)...),Vec(fill(NaN,d)...))
convert(::Type{PeriodicInterval{T}},::AnyDomain) where {T} = PeriodicInterval(nan(T),nan(T))


## Information

Base.first(d::PeriodicInterval) = d.a

Base.issubset(a::PeriodicInterval,b::PeriodicInterval) = first(a)∈b && (a.b∈b || a.b == b.b)

# we disable last since the domain is "periodic"


## Map periodic interval


tocanonical(d::PeriodicInterval{T},x) where {T} = π*(tocanonical(Segment(d),x)+1)
tocanonicalD(d::PeriodicInterval{T},x) where {T} = π*tocanonicalD(Segment(d),x)
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
