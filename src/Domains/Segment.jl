export Segment


convert(::Type{ChebyshevInterval}, ::AnyDomain) = ChebyshevInterval()
convert(::Type{ChebyshevInterval{T}}, ::AnyDomain) where T = ChebyshevInterval{T}()


## Standard interval
# T Must be a Vector space
"""
	Segment(a,b)

represents a line segment from `a` to `b`.  In the case where `a` and `b`
are real and `a < b`, then this is is equivalent to an `Interval(a,b)`.
"""
struct Segment{T} <: AbstractSegment{T}
	a::T
	b::T
	Segment{T}() where {T} = new{T}(-one(T),one(T))
	Segment{T}(a,b) where {T} = new{T}(a,b)
end


Segment() = Segment{Float64}()
Segment(a::Complex{IT1},b::Complex{IT2}) where {IT1<:Integer,IT2<:Integer} =
	Segment(ComplexF64(a),ComplexF64(b)) #convenience method
Segment(a::Integer,b::Integer) = Segment(Float64(a),Float64(b)) #convenience method
Segment(a::Complex{IT},b) where {IT<:Integer} = Segment(ComplexF64(a),b) #convenience method
Segment(a,b::Complex{IT}) where {IT<:Integer} = Segment(a,ComplexF64(b)) #convenience method
Segment(a,b) = Segment{promote_type(typeof(a),typeof(b))}(a,b)
Segment(a::Tuple,b::Tuple) = Segment(Vec(a...),Vec(b...))


convert(::Type{Domain{T}}, d::Segment) where {T<:Number} = Segment{T}(leftendpoint(d),rightendpoint(d))
convert(::Type{Segment{T}}, d::Segment) where {T<:Number} = Segment{T}(leftendpoint(d),rightendpoint(d))
convert(::Type{Segment}, d::AbstractInterval) = Segment(leftendpoint(d), rightendpoint(d))
convert(::Type{Segment{T}}, d::AbstractInterval) where T =convert(Segment{T}, convert(Segment, d))



Segment(d::AbstractInterval) = convert(Segment, d)


AnySegment(::Type{T}) where {T} = Segment{T}(NaN,NaN)
AnySegment() = AnySegment(Float64)
isambiguous(d::Segment) = all(isnan(leftendpoint(d))) && all(isnan(rightendpoint(d)))
convert(::Type{Segment{T}},::AnyDomain) where {T<:Number} = AnySegment(T)
convert(::Type{Segment},::AnyDomain) = AnySegment()


## Information
@inline leftendpoint(d::Segment) = d.a
@inline rightendpoint(d::Segment) = d.b
@inline endpoints(d::Segment) = d.a, d.b

@inline minimum(d::Segment) = min(leftendpoint(d),rightendpoint(d))
@inline maximum(d::Segment) = max(leftendpoint(d),rightendpoint(d))

isempty(d::Segment) = isapprox(leftendpoint(d), rightendpoint(d); atol=200eps(eltype(d)))

issubset(a::Segment,b::Segment) = first(a)∈b && last(a)∈b



arclength(d::AbstractInterval) = width(d)
arclength(d::Segment) = norm(complexlength(d))
complexlength(d::IntervalOrSegment) = rightendpoint(d)-leftendpoint(d)
mean(d::IntervalOrSegment) = (rightendpoint(d)+leftendpoint(d))/2
angle(d::IntervalOrSegment) = angle(complexlength(d))
sign(d::IntervalOrSegment) = sign(complexlength(d))

## Map interval
# The first definition  is the more general

mobius(S::Space,x...) = mobius(domain(S),x...)

mobius(d::ChebyshevInterval{T},x) where {T<:Number} = x
fromcanonical(d::ChebyshevInterval{T},x) where {T<:Number} = x
fromcanonicalD(d::ChebyshevInterval{T},x) where {T<:Number} = one(x)
tocanonical(d::ChebyshevInterval{T},x) where {T<:Number} = x
tocanonicalD(d::ChebyshevInterval{T},x) where {T<:Number} = one(x)

tocanonical(d::IntervalOrSegment{T},x) where {T} = 2norm(x-leftendpoint(d))/arclength(d)-1
tocanonical(d::IntervalOrSegment{T},x::Number) where {T<:Complex} = 2norm(x-leftendpoint(d))/arclength(d)-1
mobius(d::IntervalOrSegment,x) = (2x - leftendpoint(d) - rightendpoint(d))/complexlength(d)
tocanonical(d::IntervalOrSegment{T},x) where {T<:Real} = mobius(d,x)
tocanonicalD(d::IntervalOrSegment{T},x) where {T<:Real} = 2/complexlength(d)
fromcanonical(d::IntervalOrSegment{T},x) where {T<:Number} = mean(d) + complexlength(d)x/2
fromcanonical(d::IntervalOrSegment{T},x) where {T<:Vec} = mean(d) + complexlength(d)x/2
fromcanonicalD(d::IntervalOrSegment,x) = complexlength(d) / 2



==(d::Segment, m::Segment) = (isambiguous(d) && isambiguous(m)) || (leftendpoint(d) == leftendpoint(m) && rightendpoint(d) == rightendpoint(m))
function isapprox(d::Segment, m::Segment)
    tol=10E-12
    norm(leftendpoint(d)-leftendpoint(m))<tol && norm(rightendpoint(d)-rightendpoint(m))<tol
end

for op in (:(==), :isapprox)
    @eval begin
        $op(d::Segment, m::AbstractInterval) = $op(d, Segment(m))
        $op(m::AbstractInterval, d::Segment) = $op(Segment(m), d)
    end
end



## algebra

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::Segment) = Segment($op(c,leftendpoint(d)),$op(c,rightendpoint(d)))
        $op(d::Segment,c::Number) = Segment($op(leftendpoint(d),c),$op(rightendpoint(d),c))
    end
end

broadcast(::typeof(^),c::Number,d::Segment) = Segment(c^leftendpoint(d),c^rightendpoint(d))
broadcast(::typeof(^),d::Segment,c::Number) = Segment(leftendpoint(d)^c,rightendpoint(d)^c)

/(d::Segment,c::Number) = Segment(leftendpoint(d)/c,rightendpoint(d)/c)


sqrt(d::Segment)=Segment(sqrt(leftendpoint(d)),sqrt(rightendpoint(d)))

+(d1::Segment,d2::Segment)=Segment(d1.a+d2.a,d1.b+d2.b)



## intersect/union

reverseorientation(d::IntervalOrSegment) = Segment(rightendpoint(d),leftendpoint(d))

function intersect(a::Segment{T},b::Segment{V}) where {T<:Real,V<:Real}
    if first(a) > last(a)
        intersect(reverseorientation(a),b)
    elseif first(b) > last(b)
        intersect(a,reverseorientation(b))
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


function setdiff(a::Segment{T}, b::Segment{V}) where {T<:Real,V<:Real}
    # ensure a/b are well-ordered
    if first(a) > last(a)
        intersect(reverseorientation(a),b)
    elseif first(b) > last(b)
        intersect(a,reverseorientation(b))
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
