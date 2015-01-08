

opcount(::Derivative)=1
opcount(::ConstantOperator)=1
opcount(::Conversion)=1
opcount(::Multiplication)=1


opcount(M::Union(MultiplicationWrapper,SpaceOperator,DiagonalArrayOperator))=1+opcount(M.op)
opcount(A::Union(PlusOperator,TimesOperator,InterlaceOperator,DiagonalPiecewiseOperator,SumInterlaceOperator))=1+mapreduce(opcount,+,A.ops)


add_edges!(D::Derivative,nd,M,labels)=(labels[nd]=string(nd)*":"*(D.order==1?"\$D":"\$D\^"*string(D.order))*"\$")
add_edges!(A::ConstantOperator,nd,M,labels)=(labels[nd]=string(nd)*":"*string(A.c)*"I")
add_edges!(A::Conversion,nd,M,labels)=(labels[nd]=string(nd)*":"*"C")
add_edges!(A::Multiplication,nd,M,labels)=(labels[nd]=string(nd)*":"*"M")

add_edges!(A::BandedOperator,nd,M,labels)=(labels[nd]=string(nd)*":"*string(typeof(A)))


function treeadd_edges!(str,ops,node,M,labels)
    labels[node]=str

    nd=node+1
    for op in ops
        Main.Graphs.add_edge!(M,node,nd)
        add_edges!(op,nd,M,labels)
        nd+=opcount(op)
    end
end


for (WRAP,STR) in ((:SpaceOperator,:"Sw"),(:MultiplicationWrapper,:"Mw"),(:DiagonalArrayOperator,:"DiagArray"))
    @eval add_edges!(A::$WRAP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,[A.op],nd,M,labels)
end

for (OP,STR) in ((:PlusOperator,:"+"),(:TimesOperator,:"*"),(:InterlaceOperator,:"Interlace"),
    (:DiagonalPiecewiseOperator,:"DiagPiecewise"),(:SumInterlaceOperator,:"\$\\oplus\$Interlace"))
    @eval add_edges!(A::$OP,nd,M,labels)=treeadd_edges!(string(nd)*":"*$STR,A.ops,nd,M,labels)
end


function introspect(A::BandedOperator)
    require("TikzGraphs")
    m=opcount(A)

    M=Main.Graphs.simple_graph(m)
    labels=Array(ASCIIString,m)


    add_edges!(A,1,M,labels)

    Main.TikzGraphs.plot(M,labels)
end