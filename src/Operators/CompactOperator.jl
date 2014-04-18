export CompactOperator



type CompactOperator{T<:Number} <: BandedOperator{T}
    matrix::Array{T,2}
end



function matrix_addentries!(M::Array,A::ShiftArray,kr::Range1)
    for k=kr[1]:min(size(M,1),kr[end]),j=1:size(M,2)
        A[k,j-k] += M[k,j]
    end
    
    A
end


addentries!(T::CompactOperator,A::ShiftArray,kr::Range1)=matrix_addentries!(T.matrix,A,kr)

bandrange(T::CompactOperator)=1-size(T.matrix,2):size(T.matrix,2)-1
