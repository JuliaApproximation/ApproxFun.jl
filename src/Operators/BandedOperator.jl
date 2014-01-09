

export MutableBandedOperator,givensreduce!,adaptiveqr

## General routines




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
  
    b.data[k,j-k] = x
    x
end

function fastgetindex(b,data,k::Integer,j::Integer)
    data[k,j-k + b.data.colindex]
end

function fastsetindex!(b,data,x,k::Integer,j::Integer)
    data[k,j-k + b.data.colindex] = x
end

##TODO: decide adaptive resize
function givensreduce!(B::MutableBandedOperator,v::Vector,k1::Integer,k2::Integer,j1::Integer)
    data = unsafe_view(B.data.data)

    a=fastgetindex(B,data,k1,j1);b=fastgetindex(B,data,k2,j1)
    sq=sqrt(a*a + b*b)
    a=a/sq;b=b/sq
    
    v[k1],v[k2] = a*v[k1] + b*v[k2],-b*v[k1] + a*v[k2]    
    
    
    #TODO: Assuming that left rows are already zero
    
    for j = j1:indexrange(B,k1)[end]
        B1 = fastgetindex(B,data,k1,j)
        B2 = fastgetindex(B,data,k2,j)
        
        fastsetindex!(B,data, a*B1 + b*B2,k1,j)
        fastsetindex!(B,data,-b*B1 + a*B2, k2,j)
    end
    
    for j=indexrange(B,k1)[end]+1:indexrange(B,k2)[end]
        B2 = fastgetindex(B,data,k2,j)
        fastsetindex!(B,data,a*B2,k2,j)
    end
    
    #TODO: assert that the remaining of k2 are zero
    B
end

function givensreduce!(B::MutableBandedOperator,v::Vector,k1::Range1,j1::Integer)
  for k=k1[2]:k1[end]
    givensreduce!(B,v,k1[1],k,j1)
  end
end

givensreduce!(B::MutableBandedOperator,v::Vector,j::Integer)=givensreduce!(B,v,j:(j-bandrange(B)[1]),j)


function backsubstitution!(B::MutableBandedOperator,u)
  n=length(u)
  b=bandrange(B)[end]
  
  for k=n:-1:1
    for j=k+1:min(n,k+b)
      u[k]-=B[k,j]*u[j]
    end
      
    u[k] /= B[k,k]
  end
  
  u
end

function adaptiveqr(A::BandedOperator,v)
  u=[v,zeros(100)]
  B=MutableBandedOperator(A)
  
  l = length(v) + 100  
  resizedata!(B,100)
  
  
  j=1
  b=-bandrange(B)[1]
  
  tol=.001eps()
  
  
  ##TODO: check js
  while norm(u[j:j+b-1]) > tol
    if j + b == l
      u = [u,zeros(l)]      
      l *= 2
      resizedata!(B,l)
    end

    
    givensreduce!(B,u,j)
    j+=1
  end
  
  backsubstitution!(B,u[1:j-1])
end


