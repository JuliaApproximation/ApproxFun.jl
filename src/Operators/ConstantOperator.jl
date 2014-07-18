export MultiplicationOperator


type ConstantOperator{T<:Number} <: BandedOperator{T}
    c::T
end

domainspace(M::ConstantOperator)=Any
rangespace(M::ConstantOperator)=Any

bandinds(T::ConstantOperator)=0,0

addentries!(C::ConstantOperator,A::ShiftArray,kr::Range1)=laurent_addentries!([.5C.c],A,kr)




