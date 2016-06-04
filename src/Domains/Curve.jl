

export Curve

"""
`Curve` Represents a domain defined by the image of a Fun
"""


immutable IntervalCurve{S<:Space,T} <: IntervalDomain{T}
    curve::Fun{S,T}
end

immutable PeriodicCurve{S<:Space,T} <: PeriodicDomain{T}
    curve::Fun{S,T}
end

typealias Curve{S,T} Union{IntervalCurve{S,T},PeriodicCurve{S,T}}


==(a::Curve,b::Curve)=a.curve==b.curve

for TYP in (:IntervalCurve,:PeriodicCurve)
    @eval points(c::$TYP,n::Integer)=c.curve(points(domain(c.curve),n))
end

checkpoints(d::Curve) = fromcanonical(d,checkpoints(domain(d.curve)))

for op in (:(Base.first),:(Base.last),:(Base.rand))
    @eval $op(c::Curve)=c.curve($op(domain(c.curve)))
end


canonicaldomain(c::Curve)=domain(c.curve)

fromcanonical{S<:Space,T<:Number}(d::Curve{S,T},v::AbstractArray)=eltype(d)[fromcanonical(d,vk) for vk in v]
fromcanonical{S<:Space,T<:Number}(c::Curve{S,T},x) = c.curve(x)
function tocanonical(c::Curve,x)
    rts=roots(c.curve-x)
    @assert length(rts)==1
    first(rts)
end


fromcanonicalD(c::Curve,x)=differentiate(c.curve)(x)


Base.in(x,d::Curve)=in(tocanonical(d,x),canonicaldomain(d))

Base.reverse(d::Curve)=Curve(reverseorientation(d.curve))

isambiguous(d::Curve)=length(d.curve)==0 && isambiguous(domain(d.curve))
Base.convert{S,T}(::Type{IntervalCurve{S,T}},::AnyDomain)=Fun([NaN],S(AnyDomain()))
Base.convert{S,T}(::Type{PeriodicCurve{S,T}},::AnyDomain)=Fun([NaN],S(AnyDomain()))
