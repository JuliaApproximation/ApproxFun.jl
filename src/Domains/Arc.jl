export Arc

mobius(a,b,c,d,z)=(a*z+b)./(c*z+d)
mobius(a,b,c,d,z::Number)=isinf(z)?a/c:(a*z+b)./(c*z+d)
mobiusinv(a,b,c,d,z)=mobius(d,-b,-c,a,z)


mobiusD(a,b,c,d,z)=(a*d-b*c)./(d+c*z).^2
mobiusinvD(a,b,c,d,z)=mobiusD(d,-b,-c,a,z)

immutable Arc{T<:Number,V<:Real} <: IntervalDomain{Complex{Float64}}
    center::T
    radius::V
    angles::(V,V)
end


Arc(c,r,t)=Arc{typeof(c),typeof(r)}(c,r,t)
Arc(c,r,t0,t1)=Arc(c,r,(t0,t1))


isambiguous(d::Arc)=isnan(d.center) && isnan(d.radius) && isnan(d.angles[1]) && isnan(d.angles[2])
Base.convert{T<:Number,V<:Number}(::Type{Arc{T,V}},::AnyDomain)=Arc{T,V}(NaN,NaN,(NaN,NaN))
Base.convert{IT<:Arc}(::Type{IT},::AnyDomain)=Arc(NaN,NaN,(NaN,NaN))



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



