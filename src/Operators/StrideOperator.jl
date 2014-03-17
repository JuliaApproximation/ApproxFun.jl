

export StrideOperator,StrideRowOperator




type StrideOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T}
    op::B
    rowindex::Int
    colindex::Int
    stride::Int  ##TODO: rowstride,colstride
    
    function StrideOperator(o,r,c,s)
        
        new(o,r,c,s)
    end
end

StrideOperator{T<:Number}(B::BandedOperator{T},r,c,rs)=StrideOperator{T,typeof(B)}(B,r,c,rs)

bandrange(S::StrideOperator)=(min(0,S.colindex + S.stride*bandrange(S.op)[1]):max(S.colindex + S.stride*bandrange(S.op)[end],0))


divrowrange(S,r)=fld(r[1]-S.rowindex,S.stride)+1:fld(r[end]-S.rowindex-1,S.stride)+1
#divcolumnrange(S,r)=fld(r[1]+1,S.stride):fld(r[end],S.stride)


# firststriderow(S,n)=S.stride*fld(n+S.stride-2-S.rowindex,S.stride)+S.rowindex+1
# laststriderow(S,n)=firststriderow(S,n-S.stride+1)

function addentries!{T<:Number}(S::StrideOperator{T},A::ShiftArray,kr::Range1)
    r1=divrowrange(S,rowrange(A))

    A1=ShiftArray(S.op,r1)
    
    for k=r1, j=columnrange(A1)
        A[S.stride*(k-1) + S.rowindex + 1,S.stride*j + S.colindex] = A1[k,j]
    end
    
    A
end



type StrideRowOperator{T<:Number,B<:RowOperator} <: RowOperator{T}
    op::B
    rowindex::Int
    stride::Int  
end

StrideRowOperator{T<:Number}(B::RowOperator{T},r,rs)=StrideRowOperator{T,typeof(B)}(B,r,rs)


Base.getindex{T<:Number}(op::StrideRowOperator{T},kr::Range1)=[((k-1)%op.stride==op.rowindex)?op.op[fld(k-1,op.stride)+1]:zero(T) for k=kr]

