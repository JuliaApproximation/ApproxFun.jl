

export Curve



struct IntervalCurve{S<:Space,T,VT} <: SegmentDomain{T}
    curve::Fun{S,T,VT}
end

struct PeriodicCurve{S<:Space,T,VT} <: PeriodicDomain{T}
    curve::Fun{S,T,VT}
end

"""
`Curve` Represents a domain defined by the image of a Fun.  Example
usage would be

```julia
x=Fun(1..2)
Curve(exp(im*x))  # represents an arc
```
"""
const Curve{S,T} = Union{IntervalCurve{S,T},PeriodicCurve{S,T}}


==(a::Curve, b::Curve) = a.curve == b.curve
isempty(::Curve) = false

for TYP in (:IntervalCurve,:PeriodicCurve)
    @eval points(c::$TYP,n::Integer) = c.curve.(points(domain(c.curve),n))
end

checkpoints(d::Curve) = fromcanonical.(Ref(d),checkpoints(domain(d.curve)))

for op in (:(leftendpoint),:(rightendpoint),:(rand))
    @eval $op(c::Curve)=c.curve($op(domain(c.curve)))
end


canonicaldomain(c::Curve) = domain(c.curve)

fromcanonical(c::Curve{S,T},x) where {S<:Space,T<:Number} = c.curve(x)
function tocanonical(c::Curve,x)
    rts=roots(c.curve-x)
    @assert length(rts)==1
    first(rts)
end


fromcanonicalD(c::Curve,x)=differentiate(c.curve)(x)

function indomain(x,c::Curve)
    rts=roots(c.curve-x)
    if length(rts) ≠ 1
        false
    else
        in(first(rts),canonicaldomain(c))
    end
end

reverseorientation(d::IntervalCurve) = IntervalCurve(reverseorientation(d.curve))
reverseorientation(d::PeriodicCurve) = PeriodicCurve(reverseorientation(d.curve))

isambiguous(d::Curve) = ncoefficients(d.curve)==0 && isambiguous(domain(d.curve))
convert(::Type{IntervalCurve{S,T}},::AnyDomain) where {S,T}=Fun(S(AnyDomain()),[NaN])
convert(::Type{PeriodicCurve{S,T}},::AnyDomain) where {S,T}=Fun(S(AnyDomain()),[NaN])


arclength(d::Curve) = linesum(ones(d))
