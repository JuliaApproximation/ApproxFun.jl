immutable Curve{S<:FunctionSpace} <:Domain
    curve::Fun{S}
end

==(a::Curve,b::Curve)=a.curve==b.curve

points(c::Curve,n)=c.curve[points(domain(c.curve),n)]
for op in (:(Base.first),:(Base.last),:(Base.rand))
    @eval $op(c::Curve)=c.curve[$op(domain(c.curve))]
end

fromcanonical(c::Curve,x)=c.curve[fromcanonical(domain(c.curve),x)]
function tocanonical(c::Curve,x)
    rts=roots(c.curve-x)
    @assert length(rts)==1
    tocanonical(c.curve,first(rts))
end


fromcanonicalD{D<:Interval}(c::Curve{D},x)=diff(c.curve)[fromcanonical(domain(c.curve),x)].*fromcanonicalD(domain(c.curve),0.)
fromcanonicalD(c::Curve,x)=diff(c.curve)[fromcanonical(domain(c.curve),x)].*fromcanonicalD(domain(c.curve),x)
