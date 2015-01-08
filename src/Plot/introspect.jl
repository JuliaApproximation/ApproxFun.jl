
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


function introspect(A::Union(BandedOperator,Functional,FunctionSpace))
    require("TikzGraphs")
    m=treecount(A)

    M=Main.Graphs.simple_graph(m)
    labels=Array(ASCIIString,m)


    add_edges!(A,1,M,labels)

    Main.TikzGraphs.plot(M,labels)
end



## Spaces

shortname(FS)=string(typeof(FS))
shortname(CS::CosSpace)="Cos"
shortname(CS::SinSpace)="Sin"
shortname(CS::Chebyshev)="\$T\$"  #*string(first(domain(CS)))*","*string(last(domain(CS)))*"]"
shortname{λ}(CS::Ultraspherical{λ})="\$U\^{"*string(λ)*"}\$"
shortname(CS::Taylor)="\$H\^+\$"
shortname(CS::Hardy{false})="\$H\^-\$"
shortname(::SumSpace)="\$\\oplus\$"
shortname(::PiecewiseSpace)="\$\\bigcup\$"
shortname{T}(A::ArraySpace{T,1})="["*string(length(A))*"]"
shortname{T}(A::ArraySpace{T,2})="["*string(size(A,1))*"\$\\times\$"*string(size(A,2))*"]"

treecount(S::Union(SumSpace,PiecewiseSpace))=1+mapreduce(treecount,+,S.spaces)
treecount(S::ArraySpace)=1+treecount(S.space)
treecount(::FunctionSpace)=1


add_edges!(FS::FunctionSpace,nd,M,labels)=(labels[nd]=string(nd)*":"*shortname(FS))

for (OP) in (:SumSpace,:PiecewiseSpace)
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*shortname(A),A.spaces,nd,M,labels)
end

for (OP) in (:ArraySpace,)
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*shortname(A),[A.space],nd,M,labels)
end




## Operators


treecount(::BandedOperator)=1

treecount(::Derivative)=1
treecount(::ConstantOperator)=1
treecount(::Conversion)=1
treecount(::Multiplication)=1


treecount(M::Union(MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,ConversionWrapper,SpaceOperator,DiagonalArrayOperator))=1+treecount(M.op)
treecount(A::Union(PlusOperator,TimesOperator,InterlaceOperator,DiagonalPiecewiseOperator,SumInterlaceOperator))=1+mapreduce(treecount,+,A.ops)

domainrangestr(A)=shortname(domainspace(A))*"\$\\rightarrow\$"*shortname(rangespace(A))


shortname(A::BandedOperator)=string(typeof(A))
shortname(D::Derivative)=(D.order==1?"\$D":"\$D\^"*string(D.order))*"\$:"*domainrangestr(D)
shortname(A::ConstantOperator)=string(A.c)*"I"
shortname(C::Conversion)="C:"*domainrangestr(C)
shortname(A::Multiplication)=space(A.f)==domainspace(A)==rangespace(A)?"M["*shortname(space(A.f))*"]":"M["*shortname(space(A.f))*"]:"domainrangestr(A)

shortname(A::SpaceOperator)="("*domainrangestr(A)*")"

add_edges!(A::BandedOperator,nd,M,labels)=(labels[nd]=string(nd)*":"*shortname(A))





for (WRAP,STR) in ((:MultiplicationWrapper,:"(M)"),(:DerivativeWrapper,:"(D)"),(:ConversionWrapper,"(C)"),(:DiagonalArrayOperator,:"DiagArray"))
    @eval add_edges!(A::$WRAP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,[A.op],nd,M,labels)
end



@eval add_edges!(A::SpaceOperator,nd,M,labels)=treeadd_edges!(string(nd)*":"*shortname(A),[A.op],nd,M,labels)




for (OP,STR) in ((:PlusOperator,:"+"),(:TimesOperator,:"*"),(:InterlaceOperator,:"Interlace"),
    (:DiagonalPiecewiseOperator,:"DiagPiecewise"),(:SumInterlaceOperator,:"\$\\oplus\$"))
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,A.ops,nd,M,labels)
end

