


## Allows flexible row index ranges
type BandedArray{T<:Number}
    data::ShiftArray{T}#TODO: probably should have more cols than rows
    colinds::(Int,Int)   #colinds[1]:colinds[end] is the column range
end

#Note: data may start at not the first row.  This corresponds to the first rows being identically zero

BandedArray(S::ShiftArray,r::Range)=BandedArray(S,(r[1],r[end]))
BandedArray(S::ShiftArray,cs::Integer)=BandedArray(S,(max(1,rowinds(S)[1] + columninds(S)[1]),cs))
BandedArray(S::ShiftArray)=BandedArray(S,rowinds(S)[end]+columninds(S)[end])



Base.size(B::BandedArray)=(rowrange(B.data)[end],B.colinds[end])
Base.size(B::BandedArray,k::Integer)=size(B)[k]

bandinds(B::BandedArray)=columninds(B.data)
bandrange(B::BandedArray)=columnrange(B.data)
columninds(B::BandedArray)=B.colinds
columnrange(B::BandedArray)=Range1(B.colinds...)  # the range of columns of the whole matrix

## the range of columns in row k
function columninds(B::BandedArray,k::Integer)  #k is the row
    br=columninds(B.data)  #same as bandinds
    cr=B.colinds
    max(br[1] + k,cr[1]),min(br[end]+k,cr[end])
end

columnrange(B::BandedArray,k::Integer)=Range1(columninds(B,k)...)

rowinds(B::BandedArray)=rowinds(B.data)
rowrange(B::BandedArray)=rowrange(B.data)  #the row range of the whole matrix
function rowinds(B::BandedArray,j::Integer)  #j is the column
    br=bandinds(B)
    rr=rowinds(B)
    max(j-br[end],rr[1]),min(rr[end],j-br[1])
end

rowrange(B::BandedArray,k::Integer)=Range1(rowinds(B,k)...)


function Base.sparse{T<:Number}(B::BandedArray{T})
  ind = B.data.colindex
  ret = spzeros(T,size(B,1),size(B,2))
    
  for k=rowrange(B)
    for j=columnrange(B,k)
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
        for j=columnrange(B,k)
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


function bamultiply{T}(::Type{T},A::BandedArray,B::BandedArray)
    @assert columnrange(A) == rowrange(B)

    ri=A.data.rowindex
    Aci=A.data.colindex    
    ci=Aci+B.data.colindex-1        
# 
     biA=bandinds(A)
     biB=bandinds(B)

    S = BandedArray(
            ShiftArray(zeros(T,size(A.data.data,1),biA[end]-biA[1]+1+biB[end]-biB[1]),
                        ri,ci),
            B.colinds
        )
        
    bamultiply(S,A,B)
end



macro bafastget(A, k,j)
    return :(getindex($A.data.data,$k +$A.data.rowindex,$j-$k+$A.data.colindex))
end




*{T<:Number}(A::BandedArray{T},B::BandedArray{T})=bamultiply(T,A,B)
*{T<:Number,M<:Number}(A::BandedArray{T},B::BandedArray{M})=bamultiply(T == Complex{Float64} || M == Complex{Float64} ? Complex{Float64} : Float64,A,B)



function bamultiply(S::BandedArray,A::BandedArray,B::BandedArray)    
    ri=A.data.rowindex    
    ci=S.data.colindex
    
    for k=rowrange(S)
        for j=columnrange(S,k)
          cinds=rowinds(B,j)
          rinds=columninds(A,k)        
        
          for m=max(rinds[1],cinds[1]):min(rinds[end],cinds[end])
                @inbounds S.data.data[k+ri,j-k+ci] += @bafastget(A,k,m)*@bafastget(B,m,j)
          end
        end
    end
  
    S
end




function *{T<:Number,M<:Number}(A::BandedArray{T},b::Vector{M})
    typ = (T == Complex{Float64} || M == Complex{Float64}) ? Complex{Float64} : Float64

    @assert columnrange(A) == 1:length(b)
    
    bamultiply(zeros(typ,rowrange(A)[end]),A,b)
end


function bamultiply(S::Vector,A::BandedArray,b::Vector)    
    for k=rowrange(A), j=columnrange(A,k)
        @inbounds S[k] += A[k,j]*b[j]
    end
    S
end

*(a::Number,B::BandedArray)=BandedArray(a*B.data,B.colinds)
*(B::BandedArray,a::Number)=BandedArray(B.data*a,B.colinds)
.*(a::Number,B::BandedArray)=BandedArray(a.*B.data,B.colinds)
.*(B::BandedArray,a::Number)=BandedArray(B.data.*a,B.colinds)


function +(A::BandedArray,B::BandedArray)
    @assert A.colinds == B.colinds

    BandedArray(A.data + B.data,A.colinds)
end

function -(A::BandedArray,B::BandedArray)
    @assert A.colinds == B.colinds

    BandedArray(A.data - B.data,A.colinds)
end



