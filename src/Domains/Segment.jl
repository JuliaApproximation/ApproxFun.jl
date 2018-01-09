

export Interval, Segment



## Standard interval
# T Must be a Vector space
doc"""
	Segment(a,b)

represents a line segment from `a` to `b`.  In the case where `a` and `b`
are real and `a < b`, then this is is equivalent to an `Interval(a,b)`.
"""
struct Segment{T} <: IntervalDomain{T}
	a::T
	b::T
	Segment{T}() where {T} = new{T}(-one(T),one(T))
	Segment{T}(a,b) where {T} = new{T}(a,b)
end



Segment() = Segment{Float64}()
Segment(a::Complex{IT1},b::Complex{IT2}) where {IT1<:Integer,IT2<:Integer} =
	Segment(Complex128(a),Complex128(b)) #convenience method
Segment(a::Integer,b::Integer) = Segment(Float64(a),Float64(b)) #convenience method
Segment(a::Complex{IT},b) where {IT<:Integer} = Segment(Complex128(a),b) #convenience method
Segment(a,b::Complex{IT}) where {IT<:Integer} = Segment(a,Complex128(b)) #convenience method
Segment(a,b) = Segment{promote_type(typeof(a),typeof(b))}(a,b)
Segment(a::Tuple,b::Tuple) = Segment(Vec(a...),Vec(b...))

for TYP in (:Integer,:Real)
	@eval function Interval(a::$TYP,b::$TYP)
		if b ≤ a
			error("Interval(a,b) requires a < b.  Use Segment($a,$b) to represent the line segment from $a to $b.")
		end
		Segment(a,b)
	end
end

Interval(a,b) = error("Interval(a,b) only defined for real a and b.  Use Segment($a,$b) to represent the line segment from $a to $b.")
Interval() = Segment()

doc"""
	Interval(a::Real,b::Real)

represents the set `{x : a ≤ x ≤ b}`.
"""
Interval

convert(::Type{Segment{T}}, d::Segment) where {T<:Number} = Segment{T}(d.a,d.b)
convert(::Type{Segment},d::ClosedInterval) = Segment(d.left,d.right)

AnySegment(::Type{T}) where {T} = Segment{T}(NaN,NaN)
AnySegment() = AnySegment(Float64)
isambiguous(d::Segment) = all(isnan(d.a)) && all(isnan(d.b))
convert(::Type{Segment{T}},::AnyDomain) where {T<:Number} = AnySegment(T)
convert(::Type{Segment},::AnyDomain) = AnySegment()


## Information

Base.first(d::Segment) = d.a
Base.last(d::Segment) = d.b
Base.isempty(d::Segment) = isapprox(d.a,d.b;atol=200eps(eltype(d)))

Base.in(x::S, d::Segment{T}) where {S<:Real, T<:Real} = min(d.a, d.b) ≤ x ≤ max(d.a, d.b)
Base.issubset(a::Segment,b::Segment) = first(a)∈b && last(a)∈b



## Map interval
# The first definition  is the more general

mobius(S::Space,x...) = mobius(domain(S),x...)

tocanonical(d::Segment{T},x) where {T} = 2norm(x-d.a)/arclength(d)-1
tocanonical(d::Segment{T},x::Number) where {T<:Complex} = 2norm(x-d.a)/arclength(d)-1
mobius(d::Segment,x) = (d.a + d.b - 2x)/(d.a - d.b)
tocanonical(d::Segment{T},x) where {T<:Real} = mobius(d,x)
tocanonicalD(d::Segment{T},x) where {T<:Real} = 2/(d.b- d.a)
fromcanonical(d::Segment{T},x) where {T<:Number} = (d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonical(d::Segment{T},x) where {T<:Vec} = (d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonicalD(d::Segment,x) = (d.b- d.a) / 2


arclength(d::Segment) = norm(d.b - d.a)
Base.angle(d::Segment) = angle(d.b-d.a)
Base.sign(d::Segment) = sign(d.b-d.a)
complexlength(d::Segment) = d.b-d.a


==(d::Segment,m::Segment) = (isambiguous(d) && isambiguous(m)) || (d.a == m.a && d.b == m.b)
function Base.isapprox(d::Segment,m::Segment)
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


Base.sqrt(d::Segment)=Segment(sqrt(d.a),sqrt(d.b))

+(d1::Segment,d2::Segment)=Segment(d1.a+d2.a,d1.b+d2.b)



## intersect/union

Base.reverse(d::Segment)=Segment(d.b,d.a)

function Base.intersect(a::Segment{T},b::Segment{V}) where {T<:Real,V<:Real}
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


function Base.setdiff(a::Segment{T},b::Segment{V}) where {T<:Real,V<:Real}
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
Base.isless(d1::Segment{T1},d2::Segment{T2}) where {T1<:Real,T2<:Real}=d1≤first(d2)&&d1≤last(d2)
Base.isless(d1::Segment{T},x::Real) where {T<:Real}=first(d1)≤x && last(d1)≤x
Base.isless(x::Real,d1::Segment{T}) where {T<:Real}=x≤first(d1) && x≤last(d1)
