export BandedOperator,givensreduce!

type BandedOperator <: InfiniteOperator
  op::InfiniteOperator
  data::Vector{Vector{Float64}} 
  datalength::Integer
end


BandedOperator(op::InfiniteOperator)=BandedOperator(op,[],0);


#for bandrange, we save room for changed entries during Givens
bandrange(b::BandedOperator)=bandrange(b.op)[1]:(bandrange(b.op)[end]-bandrange(b.op)[1])


function indexrange(b::BandedOperator,k::Integer)
    ret = bandrange(b) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end


datalength(b::BandedOperator)=b.datalength


function Base.getindex(b::BandedOperator,kr::Range1,jr::Range1)
  ret = spzeros(length(kr),length(jr))
  
  for k = kr
      ir = bandrange(b) + k
    
    for j=max(ir[1],jr[1]):min(ir[end],jr[end])
  
      if k > datalength(b)
        ret[k-kr[1]+1,j-jr[1]+1] = b.op[k,j]
      else
        ret[k-kr[1]+1,j-jr[1]+1] = b.data[k][j - ir[1] + 1]
      end
    end
  end
  
  ret
end

getindex!(b::BandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::BandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!(b::BandedOperator,n::Integer)
  l = datalength(b)
  if n > l
    resize!(b.data,2n)
    
    for k=l+1:n
      b.data[k] = Array(Float64,length(bandrange(b)))
    end

    sh = bandrange(b)[1]    
    
    for k=l+1:n
      for j=indexrange(b,k)
        b.data[k][j-k-sh+1] = b.op[k,j]
      end
      
      b.datalength = n;
    end
  end
  
  b
end

function Base.setindex!(b::BandedOperator,x,k::Integer,j::Integer)
  resizedata!(b,k)
  
  if abs(b[k,j] - x) > 10eps()
      sh = bandrange(b)[1] 
      b.data[k][j-k-sh+1] = x
      x
  end
end

function givensreduce!(B::BandedOperator,k1::Integer,k2::Integer,j1::Integer)
  a=getindex!(B,k1,j1);b=getindex!(B,k2,j1);
  sq=sqrt(a*a + b*b);
  a=a/sq;b=b/sq;
  
  #TODO: Assuming that left rows are already zero
  
  for j = j1:indexrange(B,k1)[end]
    na = a*B[k1,j] + b*B[k2,j]
    nb = -b*B[k1,j] + a*B[k2,j]
    
    B[k1,j] = na
    B[k2,j] = nb
  end
  
  #TODO: assert that the remaining of k2 are zero
  B
end

function givensreduce!(B::BandedOperator,k1::Range1,j1::Integer)
  for k=k1[2]:k1[end]
    givensreduce!(B,k1[1],k,j1)
  end
end

givensreduce!(B::BandedOperator,j::Integer)=givensreduce!(B,j:(j-bandrange(B)[1]),j)
