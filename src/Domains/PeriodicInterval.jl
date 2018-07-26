

export PeriodicInterval

"""
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

function convert(::Type{PeriodicInterval}, d::IntervalSets.ClosedInterval)
	a,b = d.left,d.right
    @assert isfinite(a) && isfinite(b)
    PeriodicInterval(a,b)
end

Segment(d::PeriodicInterval) = Segment(leftendpoint(d),rightendpoint(d))
Interval(d::PeriodicInterval) = Interval(leftendpoint(d),rightendpoint(d))
PeriodicInterval(d::Segment) = PeriodicInterval(leftendpoint(d),rightendpoint(d))

convert(::Type{PeriodicInterval{T}}, d::PeriodicInterval) where {T<:Number} = PeriodicInterval{T}(leftendpoint(d),rightendpoint(d))

isambiguous(d::PeriodicInterval) = all(isnan(leftendpoint(d))) && all(isnan(rightendpoint(d)))
convert(::Type{PeriodicInterval{T}},::AnyDomain) where {T<:Number} = PeriodicInterval{T}(NaN,NaN)
convert(::Type{PeriodicInterval{Vec{d,T}}},::AnyDomain) where {d,T} = PeriodicInterval(Vec(fill(NaN,d)...),Vec(fill(NaN,d)...))
convert(::Type{PeriodicInterval{T}},::AnyDomain) where {T} = PeriodicInterval(nan(T),nan(T))
PeriodicInterval{T}(d) where T = convert(PeriodicInterval{T}, d)
PeriodicInterval(d) = convert(PeriodicInterval, d)


## Information
leftendpoint(d::PeriodicInterval) = d.a
rightendpoint(d::PeriodicInterval) = d.b

first(d::PeriodicInterval) = leftendpoint(d)

issubset(a::PeriodicInterval,b::PeriodicInterval) = first(a)∈b && (a.b∈b || a.b == b.b)

# we disable last since the domain is "periodic"


## Map periodic interval


tocanonical(d::PeriodicInterval{T},x) where {T} = π*(tocanonical(Segment(d),x)+1)
tocanonicalD(d::PeriodicInterval{T},x) where {T} = π*tocanonicalD(Segment(d),x)
fromcanonical(d::PeriodicInterval,θ) = fromcanonical(Segment(d),θ/π-1)
fromcanonicalD(d::PeriodicInterval,θ) = fromcanonicalD(Segment(d),θ/π-1)/π



arclength(d::PeriodicInterval) = norm(complexlength(d))
angle(d::PeriodicInterval) = angle(complexlength(d))
complexlength(d::PeriodicInterval) = rightendpoint(d) - leftendpoint(d)
reverseorientation(d::PeriodicInterval) = PeriodicInterval(rightendpoint(d), leftendpoint(d))



==(d::PeriodicInterval,m::PeriodicInterval) = leftendpoint(d) == m.a && rightendpoint(d) == m.b





## algebra

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::PeriodicInterval) = PeriodicInterval($op(c,leftendpoint(d)),$op(c,rightendpoint(d)))
        $op(d::PeriodicInterval,c::Number) = PeriodicInterval($op(leftendpoint(d),c),$op(rightendpoint(d),c))
    end
end


@eval /(d::PeriodicInterval,c::Number) = PeriodicInterval(/(leftendpoint(d),c),/(rightendpoint(d),c))

+(d1::PeriodicInterval,d2::PeriodicInterval) = PeriodicInterval(d1.a+d2.a,d1.b+d2.b)
