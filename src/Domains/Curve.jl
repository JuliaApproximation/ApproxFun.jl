##
# Represents a domain defined by the image of a Fun
###

export Curve

immutable Curve{S<:Space,T} <: UnivariateDomain{T}
    curve::Fun{S,T}
end

==(a::Curve,b::Curve)=a.curve==b.curve

points(c::Curve,n)=c.curve[points(domain(c.curve),n)]
checkpoints(d::Curve)=fromcanonical(d,checkpoints(domain(d.curve)))

for op in (:(Base.first),:(Base.last),:(Base.rand))
    @eval $op(c::Curve)=c.curve[$op(domain(c.curve))]
end


canonicaldomain(c::Curve)=domain(c.curve)


fromcanonical(c::Curve,x)=c.curve[x]
function tocanonical(c::Curve,x)
    rts=roots(c.curve-x)
    @assert length(rts)==1
    first(rts)
end


fromcanonicalD(c::Curve,x)=differentiate(c.curve)[x]


Base.in(x,d::Curve)=in(tocanonical(d,x),canonicaldomain(d))
