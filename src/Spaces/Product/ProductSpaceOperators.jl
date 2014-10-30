

Derivative(sp::PiecewiseSpace)=DerivativeWrapper(interlace(diagm(map(Derivative,sp.spaces))),1)
Derivative(sp::PiecewiseSpace,k::Integer)=DerivativeWrapper(interlace(diagm(map(s->Derivative(s,k),sp.spaces))),k)


function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    vf=vec(f)
    @assert length(vf)==length(sp)
    MultiplicationWrapper(f,interlace(diagm([Multiplication(vf[k],sp[k]) for k=1:length(vf)])))
end

function Conversion(f::PiecewiseSpace,g::PiecewiseSpace)
    @assert length(f)==length(g)
    ConversionWrapper(interlace(diagm(Operator[Conversion(f[k],g[k]) for k=1:length(f)])))
end

for op in (:dirichlet,:neumann)
    @eval $op(d::PiecewiseSpace)=interlace($op(d.spaces))
    @eval $op(d::UnionDomain)=interlace($op(d.domains))
end