export ConstantOperator, BasisFunctional


type ConstantOperator{T<:Union(Float64,Complex{Float64})} <: BandedOperator{T}
    c::T
end

ConstantOperator(c::Number)=ConstantOperator(1.0c)
ConstantOperator(L::UniformScaling)=ConstantOperator(L.λ)
IdentityOperator()=ConstantOperator(1.0)

domainspace(M::ConstantOperator)=AnySpace()
rangespace(M::ConstantOperator)=AnySpace()

bandinds(T::ConstantOperator)=0,0

addentries!(C::ConstantOperator,A::ShiftArray,kr::Range1)=laurent_addentries!([.5C.c],A,kr)

==(C1::ConstantOperator,C2::ConstantOperator)=C1.c==C2.c


## Algebra

for op in (:+,:-,:*)
    @eval ($op)(A::ConstantOperator,B::ConstantOperator)=ConstantOperator($op(A.c,B.c))
end


## Basis Functional

type BasisFunctional <: Functional{Float64}
    k::Integer
end


Base.getindex(op::BasisFunctional,k::Integer)=(k==op.k)?1.:0.
Base.getindex(op::BasisFunctional,k::Range1)=convert(Vector{Float64},k.==op.k)

