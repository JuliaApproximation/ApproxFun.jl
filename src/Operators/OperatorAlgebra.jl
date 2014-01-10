

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
    ops::Vector{BandedOperator}
end

bandrange(P::TimesOperator)=mapreduce(x->bandrange(x)[1],+,P.ops):mapreduce(x->bandrange(x)[end],+,P.ops)

function addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    kre=kr[1]:(kr[end]+mapreduce(x->bandrange(x)[end],+,P.ops[1:end-1]))

    Z = ShiftArray(zeros(length(kre),size(A,2)),A.colindex,1-kr[1])
    addentries!(P.ops[end],Z,kre)
    
    for j=length(P.ops)-1:-1:1
        krr=kr[1]:(kr[end]+mapreduce(x->bandrange(x)[end],+,P.ops[1:j-1]))    
        multiplyentries!(P.ops[j],Z,krr)
    end
    
    for k=kr,j=columnrange(A)
        A[k,j] += Z[k,j]
    end
    
    A
end

*(A::TimesOperator,B::TimesOperator)=TimesOperator([A.ops,B.ops])
*(A::TimesOperator,B::BandedOperator)=TimesOperator([A.ops,B])
*(A::BandedOperator,B::TimesOperator)=TimesOperator([A,B.ops])
*(A::BandedOperator,B::BandedOperator)=TimesOperator([A,B])