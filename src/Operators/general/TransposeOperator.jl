

export TransposeOperator




mutable struct TransposeOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
end

TransposeOperator(B::Operator{T}) where {T<:Number}=TransposeOperator{T,typeof(B)}(B)

convert(::Type{Operator{T}},A::TransposeOperator) where {T}=TransposeOperator(convert(Operator{T},A.op))

domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)

domain(P::TransposeOperator)=domain(P.op)

bandinds(P::TransposeOperator)=-bandinds(P.op)[end],-bandinds(P.op)[1]


getindex(P::TransposeOperator,k::Integer,j::Integer) =
    P.op[j,k]

function convert(::Type{BandedMatrix},S::SubOperator{T,TO}) where {T,TO<:TransposeOperator}
    kr,jr=parentindexes(S)
    BandedMatrix(view(parent(S).op,jr,kr)).'
end


Base.transpose(A::Operator)=TransposeOperator(A)
Base.ctranspose(A::Operator{T}) where {T<:Real}=TransposeOperator(A)
