

export StrideOperator,StrideRowOperator



#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == op[(k-rowindex)/rowstride,(j-colindex)/colstride]
type StrideOperator{T<:Number,B<:Operator{T}} <: BandedOperator{T}
    op::B
    rowindex::Int       
    colindex::Int       
    rowstride::Int
    colstride::Int
    
    function StrideOperator(o,r,c,rs,cs)
        @assert abs(rs) == abs(cs)
        
        new(o,r,c,rs,cs)
    end
end

StrideOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=StrideOperator{T,typeof(B)}(B,r,c,rs,cs)
StrideOperator{T<:Number}(B::Operator{T},r,c,rs)=StrideOperator{T,typeof(B)}(B,r,c,rs,rs)

function bandrange(S::StrideOperator)
    br=bandrange(S.op)
    
    st = S.colstride
    
    if S.colstride >= 0 && S.rowstride >= 0
        (st*br[1]:st*br[end])-S.rowindex+S.colindex
    elseif S.colstride >= 0
        -(S.colindex+S.rowindex-2+bandrange(S.op)[end]):S.colindex+S.rowindex-2+bandrange(S.op)[end]        
    elseif S.rowstride >= 0
        -(S.colindex+S.rowindex-2-bandrange(S.op)[1]):S.colindex+S.rowindex-2-bandrange(S.op)[1]    
    else
        error("negative negative not implemented")
    end
end


function divrowrange(S,r)
    if S.rowstride > 0 && S.colstride > 0
        div(r[1] - S.rowindex+S.rowstride-1,S.rowstride):div(r[end]-S.rowindex,S.rowstride)
    elseif S.rowstride > 0
        div(r[1] - S.rowindex+S.rowstride-1,S.rowstride):min(div(r[end]-S.rowindex,S.rowstride),S.colindex-bandrange(S.op)[1]-1)
    elseif S.rowstride < 0
        max(-S.colindex-bandrange(S.op)[end]+1,-div(r[end] -S.rowindex,-S.rowstride)):-div(r[1]-S.rowindex,-S.rowstride)
    end
end



#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == A[k,j-k]
#A[rowstride*k + rowindex,colstride*j + colindex - k] == op[k,j]

function stride_pospos_addentries!(S::StrideOperator,A::ShiftArray,kr::Range1)
    r1=divrowrange(S,kr)

    B1=BandedArray(S.op,r1)
    B=BandedArray(A)
    
    for k=r1, j=columnrange(B1.data)+k
        B[S.rowstride*k + S.rowindex,S.colstride*j + S.colindex] = B1.data[k,j-k]
    end
    
    A
end

function stride_posneg_addentries!(S::StrideOperator,A::ShiftArray,kr::Range1)
    r1=divrowrange(S,kr)
    B1=ShiftArray(S.op,r1)
    B=BandedArray(A)
    
    for k=r1, j=bandrange(S.op)
        if S.colstride*(j+k) + S.colindex > 0
            B[S.rowstride*k + S.rowindex,S.colstride*(j+k) + S.colindex] = B1[k,j]
        end
    end

    
    A
end

function addentries!(S::StrideOperator,A,kr)
    if S.rowstride > 0 && S.colstride > 0
        stride_pospos_addentries!(S,A,kr)
    elseif S.rowstride > 0
        stride_posneg_addentries!(S,A,kr)    
    elseif S.colstride > 0
        stride_posneg_addentries!(S,A,kr)            
    else
        stride_negneg_addentries!(S,A,kr)            
    end
end


domain(S::StrideOperator)=domain(S.op)


## StrideRowOperator


type StrideRowOperator{T<:Number,B<:RowOperator} <: RowOperator{T}
    op::B
    rowindex::Int
    stride::Int  
end

StrideRowOperator{T<:Number}(B::RowOperator{T},r,rs)=StrideRowOperator{T,typeof(B)}(B,r,rs)


Base.getindex{T<:Number}(op::StrideRowOperator{T},kr::Range1)=[((k-1)%op.stride==op.rowindex)?op.op[fld(k-1,op.stride)+1]:zero(T) for k=kr]



##interlace block operators

iszerooperator(A::ConstantOperator)=A.c==0.
iszerooperator(A)=false
function isboundaryrow(A,k)
    for j=1:size(A,2)
        if typeof(A[k,j]) <: RowOperator
            return true
        end
    end
        
    return false
end

function interlace{T<:Operator}(A::Array{T,2})
    for k=1:size(A,1)-2
        @assert isboundaryrow(A,k) 
    end
    
    S=Array(Operator,size(A,1)-1)
    
    for k=1:size(A,1)-2
        if iszerooperator(A[k,2])
            S[k] = StrideRowOperator(A[k,1],0,2)
        elseif iszerooperator(A[k,1])
            S[k] = StrideRowOperator(A[k,2],1,2)
        else
            S[k] = StrideRowOperator(A[k,1],0,2)+StrideRowOperator(A[k,2],1,2)
        end
    end

    A31,A32=promotespaces([A[end-1,1],A[end-1,2]])
    A41,A42=promotespaces([A[end,1],A[end,2]])

    S[end] = StrideOperator(A31,0,0,2) + StrideOperator(A32,0,1,2) + 
    StrideOperator(A41,1,-1,2) + StrideOperator(A42,1,0,2) 
    
    if(size(S,1) ==1)
        S[1]
    else
        S
    end
end

## only works for BandedShiftOperator
function interlace(A::Operator)
    S1=StrideOperator(A,2,0,2)
    S2=StrideOperator(A,1,0,-2)
    S1+S2
end


