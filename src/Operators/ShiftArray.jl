export ShiftArray,BandedArray

type ShiftArray{T<:Number}
  data::Array{T,2}
  colindex::Int
  rowindex::Int
end

ShiftArray(dat::Array,ind::Integer)=ShiftArray(dat,ind,0)
ShiftArray(ind::Integer)=ShiftArray(Array(Float64,0,0),ind,0)

Base.setindex!{T<:Number}(S::ShiftArray{T},x::T,k::Integer,j::Integer)=(S.data[k + S.rowindex, j + S.colindex] = x)

Base.getindex(S::ShiftArray,k,j) = S.data[k + S.rowindex, j + S.colindex]


Base.size(S::ShiftArray,k)=size(S.data,k)
columnrange(A::ShiftArray)=(1:size(A,2))-A.colindex
rowrange(A::ShiftArray)=(1:size(A,1))-A.rowindex

*(S::ShiftArray,x::Number)=ShiftArray(x*S.data,S.colindex,S.rowindex)

function Base.resize!(S::ShiftArray,n::Integer,m::Integer)
    olddata = S.data
    S.data = zeros(n,m)
    
    if size(olddata,1) != 0
        S.data[1:size(olddata,1),1:m] = olddata[:,1:m]
    end
    
    S
end



type BandedArray{T<:Number}
    data::ShiftArray{T}#TODO: probably should have more cols than rows
end

function BandedArray(B::BandedOperator,k::Range1)
    A=ShiftArray(zeros(length(k),length(bandrange(B))),index(B))
    BandedArray(addentries!(B,A,k))
end



Base.size(B::BandedArray)=(size(B.data.data,1),size(B.data.data,1)+columnrange(B.data)[end])
Base.size(B::BandedArray,k::Integer)=size(B)[k]

bandrange(B::BandedArray)=1-B.data.colindex:size(B.data,2) - B.data.colindex
function indexrange(B::BandedArray,k::Integer)
    ret = bandrange(B) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end

columnindexrange(B::BandedArray,j::Integer)=max(1,j-bandrange(B)[end]):min((j-bandrange(B)[1]),size(B,1))

function Base.sparse(B::BandedArray)
  ind = B.data.colindex
  ret = spzeros(size(B.data,1),size(B.data,1)+size(B.data,2)-ind)
    
  for k=1:size(B.data,1)
    for j=indexrange(B,k)
      ret[k,j]=B.data[k,j-k]
    end
  end
    
    ret
end

Base.full(B::BandedArray)=full(sparse(B))

Base.getindex(B::BandedArray,k,j)=sparse(B)[k,j]
Base.setindex!(B::BandedArray,x,k::Integer,j::Integer)=(B.data[k,j-k]=x)


function *(A::BandedArray,B::BandedArray)
    @assert size(A,2) == size(B,1)
  
    S = BandedArray(ShiftArray(zeros(size(A.data.data,1),length(bandrange(A))+length(bandrange(B))),
    A.data.colindex+B.data.colindex));
    
    
    for k=1:size(S,1)
        for j=indexrange(S,k)
          for m=1:min(indexrange(A,k)[end],columnindexrange(B,j)[end])
                S[k,j] += A[k,m]*B[m,j]
          end
        end
    end
  
    S
end

