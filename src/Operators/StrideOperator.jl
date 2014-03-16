

export StrideOperator,StrideRowOperator




type StrideOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T}
    op::B
    rowindex::Int
    stride::Int  ##TODO: rowstride,colstride
end

StrideOperator{T<:Number}(B::BandedOperator{T},r,rs)=StrideOperator{T,typeof(B)}(B,r,rs)

bandrange(S::StrideOperator)=(S.stride*bandrange(S.op)[1]:S.stride*bandrange(S.op)[end])


divrowrange(S,r)=fld(r[1]-S.rowindex,S.stride)+1:fld(r[end]-S.rowindex-1,S.stride)+1
divcolumnrange(S,r)=fld(r[1]+1,S.stride):fld(r[end],S.stride)


# firststriderow(S,n)=S.stride*fld(n+S.stride-2-S.rowindex,S.stride)+S.rowindex+1
# laststriderow(S,n)=firststriderow(S,n-S.stride+1)

function addentries!{T<:Number}(S::StrideOperator{T},A::ShiftArray,kr::Range1)
    r1=divrowrange(S,rowrange(A))
    c1=divcolumnrange(S,columnrange(A))
    A1=sazeros(T,r1,c1)
    
    addentries!(S.op,A1,r1)

    ##TODO: Does this work?
    
    for k=r1, j=c1
        A[S.stride*(k-1) + S.rowindex + 1,S.stride*j] = A1[k,j]
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

