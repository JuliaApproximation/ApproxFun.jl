export MultiplicationOperator


type ConstantOperator{T<:Union(Float64,Complex{Float64})} <: BandedOperator{T}
    c::T
end

ConstantOperator(c::Number)=ConstantOperator(1.0c)

domainspace(M::ConstantOperator)=Any
rangespace(M::ConstantOperator)=Any

bandinds(T::ConstantOperator)=0,0

addentries!(C::ConstantOperator,A::ShiftArray,kr::Range1)=laurent_addentries!([.5C.c],A,kr)

==(C1::ConstantOperator,C2::ConstantOperator)=C1.c==C2.c