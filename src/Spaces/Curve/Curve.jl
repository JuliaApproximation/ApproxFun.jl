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


immutable CurveSpace{S<:FunctionSpace,D<:FunctionSpace} <:DomainSpace{Complex{Float64}}
    space::S
    domain::Curve{D}
end

CurveSpace(c::Fun)=CurveSpace(c.space,Curve(c))
Space(c::Curve)=CurveSpace(c.curve.space,c)
transform(C::CurveSpace,vals)=transform(C.space,vals)
function evaluate{C<:CurveSpace,T}(f::Fun{C,T},x::Number)
    c=f.space
    rts=roots(domain(f).curve-x)
    @assert length(rts)==1
    evaluate(Fun(f.coefficients,c.space),first(rts))
end

