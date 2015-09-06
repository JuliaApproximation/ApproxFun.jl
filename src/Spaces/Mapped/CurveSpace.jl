## CurveSpace

function evaluate{C<:CurveSpace,T}(f::Fun{C,T},x::Number)
    c=f.space
    rts=roots(domain(f).curve-x)
    @assert length(rts)==1
    evaluate(Fun(f.coefficients,c.space),first(rts))
end


identity_fun{S}(d::CurveSpace{S})=Fun(d.domain.curve.coefficients,d)


