export ShiftArray,BandedArray

type ShiftArray{T<:Number}
  data::Array{T,2}
  rowindex::Int    
  colindex::Int
end

ShiftArray(dat::Array,ind::Integer)=ShiftArray(dat,0,ind)
ShiftArray(T::DataType,ind::Integer)=ShiftArray(Array(T,0,0),0,ind)
ShiftArray(ind::Integer)=ShiftArray(Float64,ind)


function ShiftArray{T<:Number}(B::Operator{T},k::Range1,j::Range1)
    A=ShiftArray(zeros(T,length(k),length(j)),1-k[1],1-j[1])
    addentries!(B,A,k)
end


ShiftArray(B::Operator,k::Range1)=ShiftArray(B,k,bandrange(B))


Base.setindex!{T<:Number}(S::ShiftArray{T},x::T,k::Integer,j::Integer)=(S.data[k + S.rowindex, j + S.colindex] = x)

Base.getindex(S::ShiftArray,k,j) = S.data[k + S.rowindex, j + S.colindex]


Base.size(S::ShiftArray,k)=size(S.data,k)
columnrange(A::ShiftArray)=(1:size(A,2))-A.colindex
rowrange(A::ShiftArray)=(1:size(A,1))-A.rowindex

*(S::ShiftArray,x::Number)=ShiftArray(x*S.data,S.rowindex,S.colindex)
*(x::Number,S::ShiftArray)=ShiftArray(x*S.data,S.rowindex,S.colindex)
.*(S::ShiftArray,x::Number)=ShiftArray(x*S.data,S.rowindex,S.colindex)
.*(x::Number,S::ShiftArray)=ShiftArray(x*S.data,S.rowindex,S.colindex)



function shiftarray_const_addentries!(B::ShiftArray,c::Number,A::ShiftArray,kr::Range1)    
    for k=kr,j=columnrange(B)
        A[k,j] += c*B[k,j]
    end
    
    A
end



function +{T<:Number}(A::ShiftArray{T},B::ShiftArray{T})
    @assert size(A,1) == size(B,1)
    @assert A.rowindex == B.rowindex
    
    cmin=min(columnrange(A)[1],columnrange(B)[1])
    cmax=max(columnrange(A)[end],columnrange(B)[end])    
    cind=1-cmin
    
    ret=zeros(T,size(A,1),cmax-cmin+1)
    
    ret[:,columnrange(A)+cind]=A.data
    ret[:,columnrange(B)+cind]+=B.data
    
    ShiftArray(ret,A.rowindex,cind)
end

function -{T<:Number}(A::ShiftArray{T},B::ShiftArray{T})
    @assert size(A,1) == size(B,1)
    @assert A.rowindex == B.rowindex
    
    cmin=min(columnrange(A)[1],columnrange(B)[1])
    cmax=max(columnrange(A)[end],columnrange(B)[end])    
    cind=1-cmin
    
    ret=zeros(T,size(A,1),cmax-cmin+1)
    
    ret[:,columnrange(A)+cind]=A.data
    ret[:,columnrange(B)+cind]-=B.data
    
    ShiftArray(ret,A.rowindex,cind)
end


function Base.resize!{T<:Number}(S::ShiftArray{T},n::Integer,m::Integer)
    olddata = S.data
    S.data = zeros(T,n,m)
    
    if size(olddata,1) != 0
        S.data[1:size(olddata,1),1:m] = olddata[:,1:m]
    end
    
    S
end

sazeros(T::DataType,n::Range1,m::Range1)=ShiftArray(zeros(T,length(n),length(m)),1-n[1],1-m[1])
sazeros(n::Range1,m::Range1)=sazeros(Float64,n,m)



## Allows flexible row index ranges
type BandedArray{T<:Number}
    data::ShiftArray{T}#TODO: probably should have more cols than rows
    colrange::Range1{Int}
end

#Note: data may start at not the first row.  This corresponds to the first rows being identically zero


BandedArray(S::ShiftArray,cs::Integer)=BandedArray(S,max(1,rowrange(S)[1] + columnrange(S)[1]):cs)
BandedArray(S::ShiftArray)=BandedArray(S,rowrange(S)[end]+columnrange(S)[end])
BandedArray(B::Operator,k::Range1)=BandedArray(B,k,(k[1]+bandrange(B)[1]):(k[end]+bandrange(B)[end]))
BandedArray(B::Operator,k::Range1,cs)=BandedArray(ShiftArray(B,k,bandrange(B)),cs)


Base.size(B::BandedArray)=(rowrange(B.data)[end],B.colrange[end])
Base.size(B::BandedArray,k::Integer)=size(B)[k]

bandrange(B::BandedArray)=columnrange(B.data)
function indexrange(B::BandedArray,k::Integer)
    br=bandrange(B)
    cr=columnrange(B)
    max(br[1] + k,cr[1]):min(br[end]+k,cr[end])
end

rowrange(B::BandedArray)=rowrange(B.data)
function columnindexrange(B::BandedArray,j::Integer)
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



