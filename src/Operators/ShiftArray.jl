export ShiftArray,BandedArray

type ShiftArray{T<:Number}
  data::Array{T,2}#TODO: probably should have more cols than rows
  colindex::Integer
  rowindex::Integer
end

ShiftArray(dat::Array,ind::Integer)=ShiftArray(dat,ind,0)
ShiftArray(ind::Integer)=ShiftArray(Array(Float64,0,0),ind,0)

Base.setindex!(S::ShiftArray,x,k::Integer,j::Integer)=(S.data[k + S.rowindex, j + S.colindex] = x)

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



type BandedArray
    data::ShiftArray
end

Base.size(S::ShiftArray,k)=size(S.data,k)
bandrange(B::BandedArray)=1-B.data.colindex:size(B.data,2) - B.data.colindex
function indexrange(B::BandedArray,k::Integer)
    ret = bandrange(B) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end


function SparseMatrix(B::BandedArray)
  ind = B.data.colindex
  ret = spzeros(size(B.data,1),size(B.data,1)+size(B.data,2)-ind)
    
  for k=1:size(B.data,1)
    for j=indexrange(B,k)
      ret[k,j]=B.data[k,j-k]
    end
  end
    
    ret
end

Base.full(B::BandedArray)=full(SparseMatrix(B))

Base.getindex(B::BandedArray,k,j)=full(B)[k,j]
Base.setindex!(B::BandedArray,k::Integer,j::Integer)=B.data[k,j-k]
