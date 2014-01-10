

export PlusOperator,TimesOperator

type PlusOperator <: BandedOperator
    ops::Vector{BandedOperator}
end

bandrange(P::PlusOperator)=mapreduce(op->bandrange(op)[1],min,P.ops):mapreduce(op->bandrange(op)[end],max,P.ops)
differentialorder(P::PlusOperator)=mapreduce(differentialorder,max,P.ops)

function addentries!(P::PlusOperator,A::ShiftArray,kr::Range1)
    for op in P.ops
        addentries!(op,A,kr)
    end
    
    A
end


+(A::BandedOperator,B::BandedOperator)=PlusOperator([A,B])


type TimesOperator <: BandedOperator
    A::BandedOperator
    B::BandedOperator
end

bandrange(P::TimesOperator)=(bandrange(P.A)[1]+bandrange(P.B)[1]):(bandrange(P.A)[end]+bandrange(P.B)[end])

function addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    Z = ShiftArray(zeros(length(kr)+bandrange(P.A)[end],size(A,2)),A.colindex,1-kr[1])
    addentries!(P.B,Z,kr[1]:(kr[end]+bandrange(P.A)[end]))
    multiplyentries!(P.A,Z,kr)
    
    for k=kr,j=columnrange(A)
        A[k,j] += Z[k,j]
    end
    
    A
end


*(A::BandedOperator,B::BandedOperator)=TimesOperator(A,B)