

export TransposeOperator




type TransposeOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
end

TransposeOperator{T<:Number}(B::Operator{T})=TransposeOperator{T,typeof(B)}(B)

Base.convert{T}(::Type{Operator{T}},A::TransposeOperator)=TransposeOperator(convert(Operator{T},A.op))

domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)

domain(P::TransposeOperator)=domain(P.op)

bandinds(P::TransposeOperator)=-bandinds(P.op)[end],-bandinds(P.op)[1]


getindex(P::TransposeOperator,k::Integer,j::Integer) =
    P.op[j,k]

function Base.convert{T,TO<:TransposeOperator}(::Type{BandedMatrix},S::SubOperator{T,TO})
    kr,jr=parentindexes(S)
    BandedMatrix(view(parent(S).op,jr,kr)).'
end


Base.transpose(A::Operator)=TransposeOperator(A)
Base.ctranspose{T<:Real}(A::Operator{T})=TransposeOperator(A)
