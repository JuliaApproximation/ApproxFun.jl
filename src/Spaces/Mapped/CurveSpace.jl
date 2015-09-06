## CurveSpace

function evaluate{C<:Curve,S,V,T}(f::Fun{MappedSpace{S,C,V},T},x::Number)
    c=f.space
    rts=roots(domain(f).curve-x)
    @assert length(rts)==1
    evaluate(Fun(f.coefficients,c.space),first(rts))
end


identity_fun{C<:Curve,S,V}(d::MappedSpace{S,C,V})=Fun(d.domain.curve.coefficients,d)


