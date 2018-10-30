

export PeriodicSegment

"""
	PeriodicSegment(a,b)

represents a periodic interval from `a` to `b`, that is, the point
`b` is identified with `a`.
"""
struct PeriodicSegment{T} <: PeriodicDomain{T}
    a::T
    b::T
    PeriodicSegment{T}() where {T} = new{T}(0,2convert(T,π))
    PeriodicSegment{T}(a,b) where {T} = new{T}(a,b)
end

PeriodicSegment()=PeriodicSegment{Float64}()
PeriodicSegment(a::Int,b::Int) = PeriodicSegment(Float64(a),Float64(b)) #convenience method
PeriodicSegment(a,b) = PeriodicSegment{promote_type(typeof(a),typeof(b))}(a,b)
PeriodicSegment(a::Tuple,b::Tuple) = Interval(Vec(a...),Vec(b...))

function convert(::Type{PeriodicSegment}, d::ClosedInterval)
	a,b = d.left,d.right
    @assert isfinite(a) && isfinite(b)
    PeriodicSegment(a,b)
end

Segment(d::PeriodicSegment) = Segment(leftendpoint(d),rightendpoint(d))
Interval(d::PeriodicSegment) = Interval(leftendpoint(d),rightendpoint(d))
PeriodicSegment(d::Segment) = PeriodicSegment(leftendpoint(d),rightendpoint(d))

convert(::Type{PeriodicSegment{T}}, d::PeriodicSegment) where {T<:Number} = PeriodicSegment{T}(leftendpoint(d),rightendpoint(d))

isambiguous(d::PeriodicSegment) = all(isnan(leftendpoint(d))) && all(isnan(rightendpoint(d)))
convert(::Type{PeriodicSegment{T}},::AnyDomain) where {T<:Number} = PeriodicSegment{T}(NaN,NaN)
convert(::Type{PeriodicSegment{Vec{d,T}}},::AnyDomain) where {d,T} = PeriodicSegment(Vec(fill(NaN,d)...),Vec(fill(NaN,d)...))
convert(::Type{PeriodicSegment{T}},::AnyDomain) where {T} = PeriodicSegment(nan(T),nan(T))
PeriodicSegment{T}(d) where T = convert(PeriodicSegment{T}, d)
PeriodicSegment(d) = convert(PeriodicSegment, d)


## Information
leftendpoint(d::PeriodicSegment) = d.a
rightendpoint(d::PeriodicSegment) = d.b

first(d::PeriodicSegment) = leftendpoint(d)

issubset(a::PeriodicSegment,b::PeriodicSegment) = first(a)∈b && (a.b∈b || a.b == b.b)

# we disable last since the domain is "periodic"


## Map periodic interval


tocanonical(d::PeriodicSegment{T},x) where {T} = π*(tocanonical(Segment(d),x)+1)
tocanonicalD(d::PeriodicSegment{T},x) where {T} = π*tocanonicalD(Segment(d),x)
fromcanonical(d::PeriodicSegment,θ) = fromcanonical(Segment(d),θ/π-1)
fromcanonicalD(d::PeriodicSegment,θ) = fromcanonicalD(Segment(d),θ/π-1)/π



arclength(d::PeriodicSegment) = norm(complexlength(d))
angle(d::PeriodicSegment) = angle(complexlength(d))
complexlength(d::PeriodicSegment) = rightendpoint(d) - leftendpoint(d)
reverseorientation(d::PeriodicSegment) = PeriodicSegment(rightendpoint(d), leftendpoint(d))



==(d::PeriodicSegment,m::PeriodicSegment) = leftendpoint(d) == m.a && rightendpoint(d) == m.b





## algebra

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::PeriodicSegment) = PeriodicSegment($op(c,leftendpoint(d)),$op(c,rightendpoint(d)))
        $op(d::PeriodicSegment,c::Number) = PeriodicSegment($op(leftendpoint(d),c),$op(rightendpoint(d),c))
    end
end


@eval /(d::PeriodicSegment,c::Number) = PeriodicSegment(/(leftendpoint(d),c),/(rightendpoint(d),c))

+(d1::PeriodicSegment,d2::PeriodicSegment) = PeriodicSegment(d1.a+d2.a,d1.b+d2.b)
