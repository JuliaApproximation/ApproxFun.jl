

export StrideOperator,StrideRowOperator




type StrideOperator{T<:Number,B<:Operator{T}} <: BandedOperator{T}
    op::B
    rowindex::Int
    colindex::Int
    stride::Int  ##TODO: rowstride,colstride
    
    function StrideOperator(o,r,c,s)
        
        new(o,r,c,s)
    end
end

StrideOperator{T<:Number}(B::Operator{T},r,c,rs)=StrideOperator{T,typeof(B)}(B,r,c,rs)

function bandrange(S::StrideOperator)
    if S.stride > 0
        min(0,S.colindex + S.stride*bandrange(S.op)[1]):max(S.colindex + S.stride*bandrange(S.op)[end],0)
    else
        min(S.colindex + S.stride*bandrange(S.op)[end],0):max(0,S.colindex + S.stride*bandrange(S.op)[1])
    end
end


function divrowrange(S,r)
    if S.stride > 0
        fld(r[1]-S.rowindex,S.stride)+1:fld(r[end]-S.rowindex-1,S.stride)+1
    else
        -fld(r[end]-S.rowindex-1,-S.stride)-1:-fld(r[1]-S.rowindex,-S.stride)-1
    end
end
#divcolumnrange(S,r)=fld(r[1]+1,S.stride):fld(r[end],S.stride)


# firststriderow(S,n)=S.stride*fld(n+S.stride-2-S.rowindex,S.stride)+S.rowindex+1
# laststriderow(S,n)=firststriderow(S,n-S.stride+1)

function addentries!{T<:Number}(S::StrideOperator{T},A::ShiftArray,kr::Range1)
    r1=divrowrange(S,rowrange(A))

    A1=ShiftArray(S.op,r1)
    
    for k=r1, j=columnrange(A1)
        if S.stride > 0
            A[S.stride*(k-1) + S.rowindex + 1,S.stride*j + S.colindex] = A1[k,j]
        else
            A[S.stride*(k+1) + S.rowindex + 1,S.stride*j + S.colindex] = A1[k,j]          
        end
    end
    
    A
end


domain(S::StrideOperator)=domain(S.op)


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


