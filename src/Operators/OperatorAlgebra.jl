

export PlusOperator,TimesOperator

type PlusOperator <: BandedOperator
    ops::Vector{BandedOperator}
end

bandrange(P::PlusOperator)=mapreduce(op->bandrange(op)[1],min,P.ops):mapreduce(op->bandrange(op)[end],max,P.ops)

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

bandrange(P::PlusOperator)=(bandrange(A)[1]+bandrange(B)[1]):(bandrange(A)[end]+bandrange(B)[end])

function addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
#     ShiftArray(zeros(n+2,
#     addentries!(M.B,A,1:n+2)
#     multiplyentries!(M.A,A,kr)

    A
end


+(A::BandedOperator,B::BandedOperator)=PlusOperator([A,B])