

export TransposeOperator




type TransposeOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T}
    op::B
end

TransposeOperator{T<:Number}(B::BandedOperator{T})=TransposeOperator{T,typeof(B)}(B)

Base.convert{T}(::Type{BandedOperator{T}},A::TransposeOperator)=TransposeOperator(convert(BandedOperator{T},A.op))

domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)

domain(P::TransposeOperator)=domain(P.op)

bandinds(P::TransposeOperator)=-bandinds(P.op)[end],-bandinds(P.op)[1]


getindex(P::TransposeOperator,k::Integer,j::Integer) =
    P.op[j,k]


Base.transpose(A::BandedOperator)=TransposeOperator(A)
Base.ctranspose{T<:Real}(A::BandedOperator{T})=TransposeOperator(A)
