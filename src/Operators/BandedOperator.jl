

export MutableBandedOperator,givensreduce!

## General routines

function indexrange(b::BandedOperator,k::Integer)
    ret = bandrange(b) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end


## MutableBandedOperator


type MutableBandedOperator{T<:Number,M<:BandedOperator} <: BandedOperator
  op::M
  data::ShiftArray{T}  
  datalength::Integer
end


MutableBandedOperator(op::BandedOperator)=MutableBandedOperator(op,ShiftArray(index(op)),0)


index(B::MutableBandedOperator)=index(B.op)


#for bandrange, we save room for changed entries during Givens
bandrange(b::MutableBandedOperator)=bandrange(b.op)[1]:(bandrange(b.op)[end]-bandrange(b.op)[1])





datalength(b::MutableBandedOperator)=b.datalength


function Base.getindex(b::MutableBandedOperator,kr::Range1,jr::Range1)
    ret = spzeros(length(kr),length(jr))
    
    for k = kr
        ir = bandrange(b) + k
        
        for j=max(ir[1],jr[1]):min(ir[end],jr[end])
            ret[k-kr[1]+1,j-jr[1]+1] = b[k,j]
        end
    end
    
    ret
end

function Base.getindex(b::MutableBandedOperator,k::Integer,j::Integer)  
    ir = bandrange(b) + k
    
    if k > datalength(b)
        b.op[k,j]
    else
        b.data[k,j - k]
    end
end

getindex!(b::MutableBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::MutableBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator}(b::MutableBandedOperator{T,M},n::Integer)
    l = datalength(b)
    if n > l
        resize!(b.data,2n,length(bandrange(b)))
        
        copybandedentries(b.op,b.data,l+1:n,bandrange(b))
        
        b.datalength = n
    end
    
    b
end

function Base.setindex!(b::MutableBandedOperator,x,k::Integer,j::Integer)
    resizedata!(b,k)
  
    sh = bandrange(b)[1] 
    b.data[k,j-k] = x
    x
end

function givensreduce!(B::MutableBandedOperator,k1::Integer,k2::Integer,j1::Integer)
    a=getindex!(B,k1,j1);b=getindex!(B,k2,j1);
    sq=sqrt(a*a + b*b);
    a=a/sq;b=b/sq;
    
    #TODO: Assuming that left rows are already zero
    
    for j = j1:indexrange(B,k1)[end]
        B1 = B[k1,j]
        B2 = B[k2,j]
        
        B[k1,j] = a*B1 + b*B2
        B[k2,j] = -b*B1 + a*B2
    end
    
    for j=indexrange(B,k1)[end]+1:indexrange(B,k2)[end]
        B[k2,j] = a*B[k2,j]
    end
    
    #TODO: assert that the remaining of k2 are zero
    B
end

function givensreduce!(B::MutableBandedOperator,k1::Range1,j1::Integer)
  for k=k1[2]:k1[end]
    givensreduce!(B,k1[1],k,j1)
  end
end

givensreduce!(B::MutableBandedOperator,j::Integer)=givensreduce!(B,j:(j-bandrange(B)[1]),j)
