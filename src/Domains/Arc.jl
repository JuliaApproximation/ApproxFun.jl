export Arc

mobius(a,b,c,d,z)=(a*z+b)./(c*z+d)
mobius(a,b,c,d,z::Number)=isinf(z)?a/c:(a*z+b)./(c*z+d)
mobiusinv(a,b,c,d,z)=mobius(d,-b,-c,a,z)


mobiusD(a,b,c,d,z)=(a*d-b*c)./(d+c*z).^2
mobiusinvD(a,b,c,d,z)=mobiusD(d,-b,-c,a,z)

immutable Arc{T<:Number,V<:Real} <: IntervalDomain{Complex{V}}
    center::T
    radius::V
    angles::Tuple{V,V}
end


Arc{T<:Number,V<:Real,V1<:Real,V2<:Real}(c::T,r::V,t::Tuple{V1,V2})=Arc{promote_type(T,V,V1,V2),promote_type(real(T),V,V1,V2)}(c,r,t)
Arc(c,r,t0,t1)=Arc(c,r,(t0,t1))


isambiguous(d::Arc)=isnan(d.center) && isnan(d.radius) && isnan(d.angles[1]) && isnan(d.angles[2])
Base.convert{T<:Number,V<:Number}(::Type{Arc{T,V}},::AnyDomain)=Arc{T,V}(NaN,NaN,(NaN,NaN))
Base.convert{IT<:Arc}(::Type{IT},::AnyDomain)=Arc(NaN,NaN,(NaN,NaN))

Base.length(d::Arc) = d.radius*(d.angles[2]-d.angles[1])


function mobiuspars(z0,r,t0,t1)
    c=exp(im*t0/2)+exp(im*t1/2)
    f=exp(im*t0/2)-exp(im*t1/2)
    d=exp(im*t0/2+im*t1/2)*r
    -c,c*(z0+d),f,f*(d-z0)
end

mobiuspars(a::Arc)=mobiuspars(a.center,a.radius,a.angles...)

for OP in (:mobius,:mobiusinv,:mobiusD,:mobiusinvD)
    @eval $OP(a::Arc,z)=$OP(mobiuspars(a)...,z)
end


tocanonical(a::Arc,x)=mobius(a,x)
tocanonicalD(a::Arc,x)=mobiusD(a,x)
fromcanonical(a::Arc,x)=mobiusinv(a,x)
fromcanonicalD(a::Arc,x)=mobiusinvD(a,x)


## information

Base.issubset(a::Arc,b::Arc)=a.center==b.center && a.radius==b.radius && (b.angles[1]≤a.angles[1]≤a.angles[2]≤b.angles[2] ||
                                                                             b.angles[1]≥a.angles[1]≥a.angles[2]≥b.angles[2])


# Algebra

for op in (:*,:.*)
    @eval begin
        $op(c::Real,d::Arc)=Arc(c*d.center,c*d.radius,(sign(c)*d.angles[1],sign(c)*d.angles[2]))
        $op(d::Arc,c::Real)=$op(c,d)
    end
end

for op in (:+,:-,:.+,:.-)
    @eval begin
        $op(c::Number,d::Arc)=Arc($op(c,d.center),d.radius,d.angles)
        $op(d::Arc,c::Number)=Arc($op(d.center,c),d.radius,d.angles)
    end
end

# allow exp(im*Interval(0,1)) for constructing arc
function Base.exp{CMP<:Complex}(d::Interval{CMP})
    @assert isapprox(real(d.a),0) && isapprox(real(d.b),0)
    Arc(0,1,(imag(d.a),imag(d.b)))
end
