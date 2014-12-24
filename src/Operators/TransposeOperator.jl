

export TransposeOperator




type TransposeOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T} 
    op::B
end

TransposeOperator{T<:Number}(B::BandedOperator{T})=TransposeOperator{T,typeof(B)}(B)


domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)

domain(P::TransposeOperator)=domain(P.op)

bandinds(P::TransposeOperator)=-bandinds(P.op)[end],-bandinds(P.op)[1]


function addentries!(P::TransposeOperator,A,kr::Range)
    br=bandinds(P.op)
    kr2=max(kr[1]-br[end],1):kr[end]-br[1]
    addentries!(P.op,IndexTranspose(A,kr),kr2).matrix
end


transpose(A::BandedOperator)=TransposeOperator(A)
ctranspose{T<:Real}(A::BandedOperator{T})=TransposeOperator(A)



