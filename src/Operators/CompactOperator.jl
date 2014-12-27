export CompactOperator



type CompactOperator{T<:Number} <: BandedOperator{T}
    matrix::Array{T,2}
end



function matrix_addentries!(M::Array,A,kr::Range)
    for k=kr[1]:min(size(M,1),kr[end]),j=1:size(M,2)
        A[k,j] += M[k,j]
    end
    
    A
end


addentries!(T::CompactOperator,A,kr::Range1)=matrix_addentries!(T.matrix,A,kr)

bandinds(T::CompactOperator)=(1-size(T.matrix,1),size(T.matrix,2)-1)
