

export TransposeOperator




mutable struct TransposeOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
end

TransposeOperator(B::Operator{T}) where {T<:Number}=TransposeOperator{T,typeof(B)}(B)

convert(::Type{Operator{T}},A::TransposeOperator) where {T}=TransposeOperator(convert(Operator{T},A.op))

domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)

domain(P::TransposeOperator)=domain(P.op)

bandwidths(P::TransposeOperator) = reverse(bandwidths(P.op))


getindex(P::TransposeOperator,k::Integer,j::Integer) =
    P.op[j,k]

function BandedMatrix(S::SubOperator{T,TO}) where {T,TO<:TransposeOperator}
    kr,jr=parentindices(S)
    transpose(BandedMatrix(view(parent(S).op,jr,kr)))
end


Base.transpose(A::Operator)=TransposeOperator(A)
adjoint(A::Operator{T}) where {T<:Real}=TransposeOperator(A)
