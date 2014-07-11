


## Allows flexible row index ranges
type BandedArray{T<:Number}
    data::ShiftArray{T}#TODO: probably should have more cols than rows
    colrange::Range1{Int}
end

#Note: data may start at not the first row.  This corresponds to the first rows being identically zero


BandedArray(S::ShiftArray,cs::Integer)=BandedArray(S,max(1,rowrange(S)[1] + columnrange(S)[1]):cs)
BandedArray(S::ShiftArray)=BandedArray(S,rowrange(S)[end]+columnrange(S)[end])



Base.size(B::BandedArray)=(rowrange(B.data)[end],B.colrange[end])
Base.size(B::BandedArray,k::Integer)=size(B)[k]

bandrange(B::BandedArray)=columnrange(B.data)
function indexrange(B::BandedArray,k::Integer)  #k is the row
    br=bandrange(B)
    cr=columnrange(B)
    max(br[1] + k,cr[1]):min(br[end]+k,cr[end])
end

rowrange(B::BandedArray)=rowrange(B.data)
function columnindexrange(B::BandedArray,j::Integer)  #j is the column
    br=bandrange(B)
    rr=rowrange(B)
    max(j-br[end],rr[1]):min(rr[end],j-br[1])
end
columnrange(B::BandedArray)=B.colrange


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

function Base.getindex{T<:Number}(B::BandedArray{T},k::Integer,j::Integer)
    B.data[k,j-k]::T
end

function Base.getindex{T}(B::BandedArray{T},kr::Range1,jr::Range1)
    ret = spzeros(T,length(kr),length(jr))
    for k=kr
        for j=indexrange(B,k)
            ret[k - kr[1] + 1, j - jr[1] + 1] = B[k,j]
        end
    end
    
    ret
end
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



function *{T<:Number}(A::BandedArray{T},B::BandedArray{T})
    @assert columnrange(A) == rowrange(B)
  
    S = BandedArray(
            ShiftArray(zeros(T,size(A.data.data,1),length(bandrange(A))+length(bandrange(B))-1),
                        A.data.rowindex,
                        A.data.colindex+B.data.colindex-1),
            B.colrange
        );
    
    
    for k=rowrange(S)
        for j=indexrange(S,k)
          for m=max(indexrange(A,k)[1],columnindexrange(B,j)[1]):min(indexrange(A,k)[end],columnindexrange(B,j)[end])
                S[k,j] += A[k,m]*B[m,j]
          end
        end
    end
  
    S
end


##TODO: Speed up: can't tell what type S is on compile time
function *{T<:Number,M<:Number}(A::BandedArray{T},B::BandedArray{M})
    typ = (T == Complex{Float64} || M == Complex{Float64}) ? Complex{Float64} : Float64

    @assert columnrange(A) == rowrange(B)
  
    S = BandedArray(
            ShiftArray(zeros(typ,size(A.data.data,1),length(bandrange(A))+length(bandrange(B))-1),
                        A.data.rowindex,
                        A.data.colindex+B.data.colindex-1),
            B.colrange
        );
    
    
    for k=rowrange(S)
        for j=indexrange(S,k)
          for m=max(indexrange(A,k)[1],columnindexrange(B,j)[1]):min(indexrange(A,k)[end],columnindexrange(B,j)[end])
                a=A[k,m]::T
                b=B[m,j]::M
                S[k,j] += a*b
          end
        end
    end
  
    S
end


function *{T<:Number,M<:Number}(A::BandedArray{T},b::Vector{M})
    typ = (T == Complex{Float64} || M == Complex{Float64}) ? Complex{Float64} : Float64

    @assert columnrange(A) == 1:length(b)
    
    ret = zeros(typ,rowrange(A)[end])
    
    for k=rowrange(A), j=indexrange(A,k)
        ret[k] += A[k,j]*b[j]
    end
    
  
    ret
end

*(a::Number,B::BandedArray)=BandedArray(a*B.data,B.colrange)
*(B::BandedArray,a::Number)=BandedArray(B.data*a,B.colrange)
.*(a::Number,B::BandedArray)=BandedArray(a.*B.data,B.colrange)
.*(B::BandedArray,a::Number)=BandedArray(B.data.*a,B.colrange)


function +(A::BandedArray,B::BandedArray)
    @assert A.colrange == B.colrange

    BandedArray(A.data + B.data,A.colrange)
end

function -(A::BandedArray,B::BandedArray)
    @assert A.colrange == B.colrange

    BandedArray(A.data - B.data,A.colrange)
end



