

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


function addentries!(P::TransposeOperator,A,kr::Range,::Colon)
    br=bandinds(P.op)
    # the number of rows we need increases when we
    # transpose
    kr2=max(kr[1]-br[end],1):kr[end]-br[1]

    B=subview(P.op,:,kr)
    for k=kr,j=max(1,k+bandinds(P,1)):k+bandinds(P,2)
        A[k,j]+=B[j,k]
    end
    A
end


Base.transpose(A::BandedOperator)=TransposeOperator(A)
Base.ctranspose{T<:Real}(A::BandedOperator{T})=TransposeOperator(A)
