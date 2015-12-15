
## General


function treeadd_edges!(str,ops,node,M,labels)
    labels[node]=str

    nd=node+1
    for op in ops
        Main.Graphs.add_edge!(M,node,nd)
        add_edges!(op,nd,M,labels)
        nd+=treecount(op)
    end
end


function introspect(A::Union{BandedOperator,Functional,Space})
    @eval import TikzGraphs
    m=treecount(A)

    M=Main.Graphs.simple_graph(m)
    labels=Array(ASCIIString,m)


    add_edges!(A,1,M,labels)

    Main.TikzGraphs.plot(M,labels)
end



## Spaces

texname(FS)=string(typeof(FS))


texname(::Line)="\$(-\\infty,\\infty)\$"
texname(d::Interval)="\$["*string(d.a)*","*string(d.b)*"]\$"



texname(CS::CosSpace)="Cos"
texname(CS::SinSpace)="Sin"
texname(CS::Chebyshev)="\$T\$"  #*string(first(domain(CS)))*","*string(last(domain(CS)))*"]"
texname{λ}(CS::Ultraspherical{λ})="\$U\^{"*string(λ)*"}\$"
texname(CS::Taylor)="\$H\^+\$"
texname(CS::Hardy{false})="\$H\^-\$"
texname(::SumSpace)="\$\\oplus\$"
texname(::PiecewiseSpace)="\$\\bigcup\$"
texname{T}(A::ArraySpace{T,1})="["*string(length(A))*"]"
texname{T}(A::ArraySpace{T,2})="["*string(size(A,1))*"\$\\times\$"*string(size(A,2))*"]"
texname(J::JacobiWeight)="\$(1+x)\^{"*string(J.α)*"}(1-x)\^{"*string(J.β)*"}\$"

treecount(S::Union{DirectSumSpace,PiecewiseSpace})=1+mapreduce(treecount,+,S.spaces)
treecount(S::Union{ArraySpace,JacobiWeight})=1+treecount(S.space)
treecount(::Space)=1


add_edges!(FS::Space,nd,M,labels)=(labels[nd]=string(nd)*":"*texname(FS))

for (OP) in (:DirectSumSpace,:PiecewiseSpace)
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*texname(A),A.spaces,nd,M,labels)
end

for (OP) in (:ArraySpace,:JacobiWeight)
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*texname(A),[A.space],nd,M,labels)
end




## Operators


treecount(::BandedOperator)=1

treecount(::Derivative)=1
treecount(::ConstantOperator)=1
treecount(::Conversion)=1
treecount(::Multiplication)=1


treecount(M::Union{MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,
                   ConversionWrapper,SpaceOperator,DiagonalArrayOperator,ConstantTimesOperator})=1+treecount(M.op)
treecount(A::Union{PlusOperator,TimesOperator,InterlaceOperator,
                   DiagonalInterlaceOperator,
                   KroneckerOperator})=1+mapreduce(treecount,+,A.ops)

domainrangestr(A)=texname(domainspace(A))*"\$\\rightarrow\$"*texname(rangespace(A))


texname(A::BandedOperator)=string(typeof(A))
texname(D::Derivative)=(D.order==1?"\$D":"\$D\^"*string(D.order))*"\$:"*domainrangestr(D)
texname(A::ConstantOperator)=string(A.c)*"I"
texname(C::Conversion)="C:"*domainrangestr(C)
texname(A::Multiplication)=space(A.f)==domainspace(A)==rangespace(A)?"M["*texname(space(A.f))*"]":"M["*texname(space(A.f))*"]:"domainrangestr(A)


texname(D::DerivativeWrapper)=(D.order==1?"\$(D":"\$(D\^"*string(D.order))*")\$"
texname(A::SpaceOperator)="("*domainrangestr(A)*")"
texname(A::ConstantTimesOperator)=string(A.c)

add_edges!(A::BandedOperator,nd,M,labels)=(labels[nd]=string(nd)*":"*texname(A))





for (WRAP,STR) in ((:MultiplicationWrapper,:"(M)"),(:ConversionWrapper,"(C)"),(:DiagonalArrayOperator,:"DiagArray"))
    @eval add_edges!(A::$WRAP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,[A.op],nd,M,labels)
end



@eval add_edges!(A::Union{SpaceOperator,DerivativeWrapper},nd,M,labels)=treeadd_edges!(string(nd)*":"*texname(A),[A.op],nd,M,labels)




for (OP,STR) in ((:PlusOperator,:"+"),(:TimesOperator,:"*"),(:InterlaceOperator,:"Interlace"),
    (:DiagonalInterlaceOperator,:"\$\\oplus\$"),
                 (:KroneckerOperator,:"\$\\otimes\$"))
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,A.ops,nd,M,labels)
end

add_edges!(A::ConstantTimesOperator,nd,M,labels)=treeadd_edges!(string(nd)*":"*string(A.c),[A.op],nd,M,labels)
