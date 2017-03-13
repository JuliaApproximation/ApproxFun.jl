

export Curve



immutable IntervalCurve{S<:Space,T} <: IntervalDomain{T}
    curve::Fun{S,T}
end

immutable PeriodicCurve{S<:Space,T} <: PeriodicDomain{T}
    curve::Fun{S,T}
end

doc"""
`Curve` Represents a domain defined by the image of a Fun.  Example
usage would be

```julia
x=Fun(1..2)
Curve(exp(im*x))  # represents an arc
```
"""
@compat const Curve{S,T} = Union{IntervalCurve{S,T},PeriodicCurve{S,T}}


==(a::Curve,b::Curve)=a.curve==b.curve

for TYP in (:IntervalCurve,:PeriodicCurve)
    @eval points(c::$TYP,n::Integer)=c.curve(points(domain(c.curve),n))
end

checkpoints(d::Curve) = fromcanonical(d,checkpoints(domain(d.curve)))

for op in (:(Base.first),:(Base.last),:(Base.rand))
    @eval $op(c::Curve)=c.curve($op(domain(c.curve)))
end


canonicaldomain(c::Curve)=domain(c.curve)

fromcanonical{S<:Space,T<:Number}(c::Curve{S,T},x) = c.curve(x)
function tocanonical(c::Curve,x)
    rts=roots(c.curve-x)
    @assert length(rts)==1
    first(rts)
end


fromcanonicalD(c::Curve,x)=differentiate(c.curve)(x)


function Base.in(x,c::Curve)
    rts=roots(c.curve-x)
    if length(rts) ≠ 1
        false
    else
        in(first(rts),canonicaldomain(c))
    end
end

Base.reverse(d::Curve) = Curve(reverseorientation(d.curve))

isambiguous(d::Curve) = ncoefficients(d.curve)==0 && isambiguous(domain(d.curve))
Base.convert{S,T}(::Type{IntervalCurve{S,T}},::AnyDomain)=Fun(S(AnyDomain()),[NaN])
Base.convert{S,T}(::Type{PeriodicCurve{S,T}},::AnyDomain)=Fun(S(AnyDomain()),[NaN])
