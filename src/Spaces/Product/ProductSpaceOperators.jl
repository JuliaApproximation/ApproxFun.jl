

Derivative(sp::PiecewiseSpace)=DerivativeWrapper(DiagonalInterlaceOperator(map(Derivative,sp.spaces)),1)
Derivative(sp::PiecewiseSpace,k::Integer)=DerivativeWrapper(DiagonalInterlaceOperator(map(s->Derivative(s,k),sp.spaces)),k)


function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    vf=vec(f)
    @assert length(vf)==length(sp)
    MultiplicationWrapper(f,DiagonalInterlaceOperator([Multiplication(vf[k],sp[k]) for k=1:length(vf)]))
end

function Conversion(f::PiecewiseSpace,g::PiecewiseSpace)
    @assert length(f)==length(g)
    ConversionWrapper(DiagonalInterlaceOperator(Operator[Conversion(f[k],g[k]) for k=1:length(f)]))
end

for op in (:dirichlet,:neumann)
    @eval $op(d::PiecewiseSpace)=interlace($op(d.spaces))
    @eval $op(d::UnionDomain)=interlace($op(d.domains))
end



## diag provides a way to convert between DiagonalInterlaceOperator and bacn

Base.diag(A::DiagonalInterlaceOperator)=A.ops
Base.diag(A::PlusOperator)=mapreduce(diag,+,A.ops)
Base.diag(A::TimesOperator)=mapreduce(diag,.*,A.ops)

for TYP in (:DerivativeWrapper,:ConversionWrapper)
    @eval Base.diag{DT<:DiagonalInterlaceOperator}(A::($TYP{DT}))=A.op.ops
end

Base.diag{FT<:PiecewiseSpace,OT<:DiagonalInterlaceOperator}(A::MultiplicationWrapper{FT,OT})=A.op.ops