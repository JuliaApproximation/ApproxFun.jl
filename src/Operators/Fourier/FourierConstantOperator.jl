export FourierConstantOperator


type FourierConstantOperator{T<:Number} <: BandedShiftOperator{T}
    c::T
end

bandrange(T::FourierConstantOperator)=0:0

function addentries!(C::FourierConstantOperator,A::ShiftArray,kr::Range1)
    for k=kr,j=columnrange(A)
        A[k,j] *= C.c
    end
    
    A
end

