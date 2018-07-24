

export Segment



## Standard interval
# T Must be a Vector space
"""
	Segment(a,b)

represents a line segment from `a` to `b`.  In the case where `a` and `b`
are real and `a < b`, then this is is equivalent to an `Interval(a,b)`.
"""
struct Segment{T} <: SegmentDomain{T}
	a::T
	b::T
	Segment{T}() where {T} = new{T}(-one(T),one(T))
	Segment{T}(a,b) where {T} = new{T}(a,b)
end


const IntervalOrSegment{T} = Union{Domains.AbstractInterval{T}, Segment{T}}



Segment() = Segment{Float64}()
Segment(a::Complex{IT1},b::Complex{IT2}) where {IT1<:Integer,IT2<:Integer} =
	Segment(ComplexF64(a),ComplexF64(b)) #convenience method
Segment(a::Integer,b::Integer) = Segment(Float64(a),Float64(b)) #convenience method
Segment(a::Complex{IT},b) where {IT<:Integer} = Segment(ComplexF64(a),b) #convenience method
Segment(a,b::Complex{IT}) where {IT<:Integer} = Segment(a,ComplexF64(b)) #convenience method
Segment(a,b) = Segment{promote_type(typeof(a),typeof(b))}(a,b)
Segment(a::Tuple,b::Tuple) = Segment(Vec(a...),Vec(b...))


convert(::Type{Domain{T}}, d::Segment) where {T<:Number} = Segment{T}(d.a,d.b)
convert(::Type{Segment{T}}, d::Segment) where {T<:Number} = Segment{T}(d.a,d.b)
convert(::Type{Segment},d::IntervalSets.ClosedInterval) = Segment(d.left,d.right)


AnySegment(::Type{T}) where {T} = Segment{T}(NaN,NaN)
AnySegment() = AnySegment(Float64)
isambiguous(d::Segment) = all(isnan(d.a)) && all(isnan(d.b))
convert(::Type{Segment{T}},::AnyDomain) where {T<:Number} = AnySegment(T)
convert(::Type{Segment},::AnyDomain) = AnySegment()


## Information
@inline Base.first(d::Segment) = d.a
@inline Base.last(d::Segment) = d.b

@inline Base.minimum(d::Segment) = min(d.a,d.b)
@inline Base.maximum(d::Segment) = max(d.a,d.b)

Base.isempty(d::Segment) = isapprox(d.a,d.b;atol=200eps(eltype(d)))

issubset(a::Segment,b::Segment) = first(a)∈b && last(a)∈b



## Map interval
# The first definition  is the more general

mobius(S::Space,x...) = mobius(domain(S),x...)

mobius(d::ChebyshevInterval{T},x) where {T<:Number} = x
fromcanonical(d::ChebyshevInterval{T},x) where {T<:Number} = x
fromcanonicalD(d::ChebyshevInterval{T},x) where {T<:Number} = one(x)
tocanonical(d::ChebyshevInterval{T},x) where {T<:Number} = x
tocanonicalD(d::ChebyshevInterval{T},x) where {T<:Number} = one(x)

tocanonical(d::IntervalOrSegment{T},x) where {T} = 2norm(x-d.a)/arclength(d)-1
tocanonical(d::IntervalOrSegment{T},x::Number) where {T<:Complex} = 2norm(x-d.a)/arclength(d)-1
mobius(d::IntervalOrSegment,x) = (d.a + d.b - 2x)/(d.a - d.b)
tocanonical(d::IntervalOrSegment{T},x) where {T<:Real} = mobius(d,x)
tocanonicalD(d::IntervalOrSegment{T},x) where {T<:Real} = 2/(d.b- d.a)
fromcanonical(d::IntervalOrSegment{T},x) where {T<:Number} = (d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonical(d::IntervalOrSegment{T},x) where {T<:Vec} = (d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonicalD(d::IntervalOrSegment,x) = (d.b- d.a) / 2


arclength(d::Domains.AbstractInterval) = norm(maximum(d) - minimum(d))
arclength(d::Segment) = norm(d.b - d.a)
Base.angle(d::IntervalOrSegment) = angle(last(d)-first(d))
Base.sign(d::IntervalOrSegment) = sign(last(d)-first(d))
complexlength(d::IntervalOrSegment) = last(d)-first(d)


==(d::Segment,m::Segment) = (isambiguous(d) && isambiguous(m)) || (d.a == m.a && d.b == m.b)
function isapprox(d::Segment,m::Segment)
    tol=10E-12
    norm(d.a-m.a)<tol&&norm(d.b-m.b)<tol
end



## algebra

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::Segment) = Segment($op(c,d.a),$op(c,d.b))
        $op(d::Segment,c::Number) = Segment($op(d.a,c),$op(d.b,c))
    end
end

broadcast(::typeof(^),c::Number,d::Segment) = Segment(c^d.a,c^d.b)
broadcast(::typeof(^),d::Segment,c::Number) = Segment(d.a^c,d.b^c)

/(d::Segment,c::Number) = Segment(d.a/c,d.b/c)


sqrt(d::Segment)=Segment(sqrt(d.a),sqrt(d.b))

+(d1::Segment,d2::Segment)=Segment(d1.a+d2.a,d1.b+d2.b)



## intersect/union

reverse(d::Segment)=Segment(d.b,d.a)

function intersect(a::Segment{T},b::Segment{V}) where {T<:Real,V<:Real}
    if first(a) > last(a)
        intersect(reverse(a),b)
    elseif first(b) > last(b)
        intersect(a,reverse(b))
    elseif first(a) > first(b)
        intersect(b,a)
    elseif last(a) <= first(b)
        EmptyDomain()
    elseif last(a)>=last(b)
        b
    elseif isapprox(first(b),last(a);atol=100eps(promote_type(T,V))/max(arclength(a),arclength(b)))
        EmptyDomain()
    else
        Segment(first(b),last(a))
    end
end


function setdiff(a::Segment{T},b::Segment{V}) where {T<:Real,V<:Real}
    # ensure a/b are well-ordered
    if first(a) > last(a)
        intersect(reverse(a),b)
    elseif first(b) > last(b)
        intersect(a,reverse(b))
    elseif first(a)< first(b)
        if last(a) <= first(b)
            a
        else # first(a) ≤ first(b) ≤last(a)
            #TODO: setdiff in the middle
            if last(a) <= last(b)
            	Segment(first(a),first(b))
			else  # first(a) ≤ first(b) ≤ last(b) ≤last(a)
				Segment(first(a),first(b))∪Segment(last(b),last(a))
			end
        end
    else #first(a)>= first(b)
        if first(a)>=last(b)
            a
        elseif last(a) <= last(b)
            EmptyDomain()
        else #first(b) < first(a) < last(b) < last(a)
            Segment(last(b),last(a))
        end
    end
end






## sort
isless(d1::Segment{T1},d2::Segment{T2}) where {T1<:Real,T2<:Real}=d1≤first(d2)&&d1≤last(d2)
isless(d1::Segment{T},x::Real) where {T<:Real}=first(d1)≤x && last(d1)≤x
isless(x::Real,d1::Segment{T}) where {T<:Real}=x≤first(d1) && x≤last(d1)
