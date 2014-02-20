export ShiftArray,BandedArray

type ShiftArray{T<:Number}
  data::Array{T,2}
  colindex::Int
  rowindex::Int  ##TODO: Change to row, col
end

ShiftArray(dat::Array,ind::Integer)=ShiftArray(dat,ind,0)
ShiftArray(T::DataType,ind::Integer)=ShiftArray(Array(T,0,0),ind,0)
ShiftArray(ind::Integer)=ShiftArray(Float64,ind)


function ShiftArray{T<:Number}(B::BandedOperator{T},k::Range1,j::Range1)
    A=ShiftArray(zeros(T,length(k),length(j)),1-j[1],1-k[1])
    addentries!(B,A,k)
end



Base.setindex!{T<:Number}(S::ShiftArray{T},x::T,k::Integer,j::Integer)=(S.data[k + S.rowindex, j + S.colindex] = x)

Base.getindex(S::ShiftArray,k,j) = S.data[k + S.rowindex, j + S.colindex]


Base.size(S::ShiftArray,k)=size(S.data,k)
columnrange(A::ShiftArray)=(1:size(A,2))-A.colindex
rowrange(A::ShiftArray)=(1:size(A,1))-A.rowindex

*(S::ShiftArray,x::Number)=ShiftArray(x*S.data,S.colindex,S.rowindex)
*(x::Number,S::ShiftArray)=ShiftArray(x*S.data,S.colindex,S.rowindex)
.*(S::ShiftArray,x::Number)=ShiftArray(x*S.data,S.colindex,S.rowindex)
.*(x::Number,S::ShiftArray)=ShiftArray(x*S.data,S.colindex,S.rowindex)


function +{T<:Number}(A::ShiftArray{T},B::ShiftArray{T})
    @assert size(A,1) == size(B,1)
    @assert A.rowindex == B.rowindex
    
    cmin=min(columnrange(A)[1],columnrange(A)[1])
    cmax=min(columnrange(A)[end],columnrange(A)[end])    
    cind=1-cmin
    
    ret=zeros(T,size(A,1),cmax-cmin+1)
    
    ret[:,columnrange(A)+cind]=A.data
    ret[:,columnrange(B)+cind]+=B.data
    
    ShiftArray(ret,cind,A.rowindex)
end


function Base.resize!{T<:Number}(S::ShiftArray{T},n::Integer,m::Integer)
    olddata = S.data
    S.data = zeros(T,n,m)
    
    if size(olddata,1) != 0
        S.data[1:size(olddata,1),1:m] = olddata[:,1:m]
    end
    
    S
end



## Allows flexible row index ranges
type BandedArray{T<:Number}
    data::ShiftArray{T}#TODO: probably should have more cols than rows
end

BandedArray(B::BandedOperator,k::Range1)=BandedArray(ShiftArray(B,k,bandrange(B)))



Base.size(B::BandedArray)=(rowrange(B.data)[end],indexrange(B,rowrange(B.data)[end])[end])
Base.size(B::BandedArray,k::Integer)=size(B)[k]

bandrange(B::BandedArray)=columnrange(B.data)
function indexrange(B::BandedArray,k::Integer)
    ret = bandrange(B) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end

rowrange(B::BandedArray)=rowrange(B.data)
columnrange(B::BandedArray)=indexrange(B,rowrange(B)[1])[1]:indexrange(B,rowrange(B)[end])[end]
columnindexrange(B::BandedArray,j::Integer)=max(1,j-bandrange(B)[end]):min((j-bandrange(B)[1]),size(B,1))

function Base.sparse{T<:Number}(B::BandedArray{T})
  ind = B.data.colindex
  ret = spzeros(T,size(B,1),size(B,2))
    
  for k=rowrange(B)
    for j=indexrange(B,k)
      ret[k,j]=B.data[k,j-k]
    end
  end
    
    ret
end

Base.full(B::BandedArray)=full(sparse(B))

function Base.getindex(B::BandedArray,k::Integer,j::Integer)
    B.data[k,j-k]
end

Base.getindex(B::BandedArray,k::Range1,j::Range1)=sparse(B)[k,j]  ##TODO: Very slow
Base.setindex!(B::BandedArray,x,k::Integer,j::Integer)=(B.data[k,j-k]=x)

# function multiplyentries!(A::BandedArray,B::BandedArray)
#     for k=rowrange(A)
#         for j=indexrange(B,k)
#           
#         
#           for m=max(indexrange(A,k)[1],columnindexrange(B,j)[1]):min(indexrange(A,k)[end],columnindexrange(B,j)[end])
#                 S[k,j] += A[k,m]*B[m,j]
#           end
#         end
#     end 
# end


function *{T<:Number,M<:Number}(A::BandedArray{T},B::BandedArray{M})
    typ = (T == Complex{Float64} || M == Complex{Float64}) ? Complex{Float64} : Float64

    @assert columnrange(A) == rowrange(B)
  
    S = BandedArray(
        ShiftArray(zeros(typ,size(A.data.data,1),length(bandrange(A))+length(bandrange(B))-1),
        A.data.colindex+B.data.colindex-1,
        A.data.rowindex));
    
    
    for k=rowrange(S)
        for j=indexrange(S,k)
          for m=max(indexrange(A,k)[1],columnindexrange(B,j)[1]):min(indexrange(A,k)[end],columnindexrange(B,j)[end])
                S[k,j] += A[k,m]*B[m,j]
          end
        end
    end
  
    S
end




