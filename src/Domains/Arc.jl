export Arc

mobius(a,b,c,d,z) = (a*z+b)/(c*z+d)
mobius(a,b,c,d,z::Number) = isinf(z) ? a/c : (a*z+b)/(c*z+d)
mobiusinv(a,b,c,d,z) = mobius(d,-b,-c,a,z)


mobiusD(a,b,c,d,z) = (a*d-b*c)/(d+c*z)^2
mobiusinvD(a,b,c,d,z) = mobiusD(d,-b,-c,a,z)
"""
    Arc(c,r,(θ₁,θ₂))

represents the arc centred at `c` with radius `r` from angle `θ₁` to `θ₂`.
"""
struct Arc{T,V<:Real,TT} <: SegmentDomain{TT}
    center::T
    radius::V
    angles::Tuple{V,V}
    Arc{T,V,TT}(c,r,a) where {T,V,TT} = new{T,V,TT}(convert(T,c),convert(V,r),Tuple{V,V}(a))
end


Arc(c::T,r::V,t::Tuple{V1,V2}) where {T<:Number,V<:Real,V1<:Real,V2<:Real} =
    Arc{promote_type(T,V,V1,V2),
        promote_type(real(T),V,V1,V2),
        Complex{promote_type(real(T),V,V1,V2)}}(c,r,t)
Arc(c::Vec{2,T},r::V,t::Tuple{V,V}) where {T<:Number,V<:Real} =
    Arc{Vec{2,promote_type(T,V)},
        promote_type(real(T),V),
        Vec{2,promote_type(real(T),V)}}(c,r,t)
Arc(c::Vec{2,T},r::V,t::Tuple{V1,V2}) where {T<:Number,V<:Real,V1<:Real,V2<:Real} =
    Arc{Vec{2,promote_type(T,V,V1,V2)},
        promote_type(real(T),V,V1,V2),
        Vec{2,promote_type(T,V,V1,V2)}}(c,r,t)

Arc(c::Tuple,r,t) = Arc(Vec(c...),r,t)
Arc(c,r,t0,t1) = Arc(c,r,(t0,t1))


complex(a::Arc{V}) where {V<:Vec} = Arc(complex(a.center...),a.radius,a.angles)

isambiguous(d::Arc) =
    isnan(d.center) && isnan(d.radius) && isnan(d.angles[1]) && isnan(d.angles[2])
convert(::Type{Arc{T,V}},::AnyDomain) where {T<:Number,V<:Real} =
    Arc{T,V}(NaN,NaN,(NaN,NaN))
convert(::Type{IT},::AnyDomain) where {IT<:Arc} =
    Arc(NaN,NaN,(NaN,NaN))

isempty(d::Arc) = false

reverseorientation(a::Arc) = Arc(a.center,a.radius,reverse(a.angles))

arclength(d::Arc) = d.radius*(d.angles[2]-d.angles[1])


function mobiuspars(z0,r,t0,t1)
    c=exp(im*t0/2)+exp(im*t1/2)
    f=exp(im*t0/2)-exp(im*t1/2)
    d=exp(im*t0/2+im*t1/2)*r
    -c,c*(z0+d),f,f*(d-z0)
end

mobiuspars(a::Arc) = mobiuspars(a.center,a.radius,a.angles...)

for OP in (:mobius,:mobiusinv,:mobiusD,:mobiusinvD)
    @eval $OP(a::Arc,z) = $OP(mobiuspars(a)...,z)
end


tocanonical(a::Arc{T},x) where {T<:Number} = real(mobius(a,x))
tocanonicalD(a::Arc{T},x) where {T<:Number} = mobiusD(a,x)
fromcanonical(a::Arc{T,V,TT},x) where {T<:Number,V<:Real,TT<:Complex} =
    mobiusinv(a,x)
fromcanonicalD(a::Arc{T},x) where {T<:Number} = mobiusinvD(a,x)

tocanonical(a::Arc{V},x::Vec) where {V<:Vec} =
    tocanonical(complex(a),complex(x...))
fromcanonical(a::Arc{V},x::Number) where {V<:Vec} =
    Vec(reim(fromcanonical(complex(a),x))...)
fromcanonicalD(a::Arc{V},x::Number) where {V<:Vec} =
    Vec(reim(fromcanonicalD(complex(a),x))...)

## information

issubset(a::Arc,b::Arc) =
    a.center==b.center && a.radius==b.radius &&
                    (b.angles[1]≤a.angles[1]≤a.angles[2]≤b.angles[2] ||
                            b.angles[1]≥a.angles[1]≥a.angles[2]≥b.angles[2])


# Algebra

*(c::Real,d::Arc) =
    Arc(c*d.center,c*d.radius,(sign(c)*d.angles[1],sign(c)*d.angles[2]))
*(d::Arc,c::Real) = *(c,d)

for op in (:+,:-)
    @eval begin
        $op(c::Number,d::Arc) = Arc($op(c,d.center),d.radius,d.angles)
        $op(d::Arc,c::Number) = Arc($op(d.center,c),d.radius,d.angles)
    end
end


# allow exp(im*Segment(0,1)) for constructing arc
function exp(d::IntervalOrSegment{<:Complex})
    @assert isapprox(real(leftendpoint(d)),0) && isapprox(real(rightendpoint(d)),0)
    Arc(0,1,(imag(leftendpoint(d)),imag(rightendpoint(d))))
end
