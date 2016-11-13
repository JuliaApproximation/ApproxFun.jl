export Arc

mobius(a,b,c,d,z) = (a*z+b)./(c*z+d)
mobius(a,b,c,d,z::Number) = isinf(z)?a/c:(a*z+b)./(c*z+d)
mobiusinv(a,b,c,d,z) = mobius(d,-b,-c,a,z)


mobiusD(a,b,c,d,z) = (a*d-b*c)./(d+c*z).^2
mobiusinvD(a,b,c,d,z) = mobiusD(d,-b,-c,a,z)

immutable Arc{T,V<:Real,TT} <: IntervalDomain{TT}
    center::T
    radius::V
    angles::Tuple{V,V}
    Arc(c,r,a) = new(T(c),V(r),Tuple{V,V}(a))
end


Arc{T<:Number,V<:Real,V1<:Real,V2<:Real}(c::T,r::V,t::Tuple{V1,V2}) =
    Arc{promote_type(T,V,V1,V2),
        promote_type(real(T),V,V1,V2),
        Complex{promote_type(real(T),V,V1,V2)}}(c,r,t)
Arc{T<:Number,V<:Real}(c::Vec{2,T},r::V,t::Tuple{V,V}) =
    Arc{Vec{2,promote_type(T,V)},
        promote_type(real(T),V),
        Vec{2,promote_type(real(T),V)}}(c,r,t)
Arc{T<:Number,V<:Real,V1<:Real,V2<:Real}(c::Vec{2,T},r::V,t::Tuple{V1,V2}) =
    Arc{Vec{2,promote_type(T,V,V1,V2)},
        promote_type(real(T),V,V1,V2),
        Vec{2,promote_type(T,V,V1,V2)}}(c,r,t)

Arc(c::Tuple,r,t) = Arc(Vec(c...),r,t)
Arc(c,r,t0,t1) = Arc(c,r,(t0,t1))


Base.complex{V<:Vec}(a::Arc{V}) = Arc(complex(a.center...),a.radius,a.angles)

isambiguous(d::Arc) =
    isnan(d.center) && isnan(d.radius) && isnan(d.angles[1]) && isnan(d.angles[2])
Base.convert{T<:Number,V<:Real}(::Type{Arc{T,V}},::AnyDomain) =
    Arc{T,V}(NaN,NaN,(NaN,NaN))
Base.convert{IT<:Arc}(::Type{IT},::AnyDomain) =
    Arc(NaN,NaN,(NaN,NaN))

Base.reverse(a::Arc) = Arc(a.center,a.radius,reverse(a.angles))

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


tocanonical{T<:Number}(a::Arc{T},x) = real(mobius(a,x))
tocanonicalD{T<:Number}(a::Arc{T},x) = mobiusD(a,x)
fromcanonical{T<:Number,V<:Real,TT<:Complex}(d::Arc{T,V,TT},v::AbstractArray) =
    eltype(d)[fromcanonical(d,vk) for vk in v]
fromcanonical{T<:Number,V<:Real,TT<:Complex}(a::Arc{T,V,TT},x) =
    mobiusinv(a,x)
fromcanonicalD{T<:Number}(a::Arc{T},x) = mobiusinvD(a,x)

tocanonical{V<:Vec}(a::Arc{V},x::Vec) =
    tocanonical(complex(a),complex(x...))
fromcanonical{V<:Vec}(a::Arc{V},x::Number) =
    Vec(reim(fromcanonical(complex(a),x))...)
fromcanonicalD{V<:Vec}(a::Arc{V},x::Number) =
    Vec(reim(fromcanonicalD(complex(a),x))...)

## information

Base.issubset(a::Arc,b::Arc) =
    a.center==b.center && a.radius==b.radius &&
                    (b.angles[1]≤a.angles[1]≤a.angles[2]≤b.angles[2] ||
                            b.angles[1]≥a.angles[1]≥a.angles[2]≥b.angles[2])


# Algebra

for op in (:*,:.*)
    @eval begin
        $op(c::Real,d::Arc) =
            Arc(c*d.center,c*d.radius,(sign(c)*d.angles[1],sign(c)*d.angles[2]))
        $op(d::Arc,c::Real) = $op(c,d)
    end
end

for op in (:+,:-,:.+,:.-)
    @eval begin
        $op(c::Number,d::Arc) = Arc($op(c,d.center),d.radius,d.angles)
        $op(d::Arc,c::Number) = Arc($op(d.center,c),d.radius,d.angles)
    end
end

# allow exp(im*Interval(0,1)) for constructing arc
function Base.exp{CMP<:Complex}(d::Interval{CMP})
    @assert isapprox(real(d.a),0) && isapprox(real(d.b),0)
    Arc(0,1,(imag(d.a),imag(d.b)))
end
