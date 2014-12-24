

## 
# Represent a banded matrix
# [ a_11 a_12 
#   a_21 a_22 a_23
#   a_31 a_32 a_33 a_34
#        a_42 a_43 a_44  ]
# ordering the data like
#       [ *     *       a_31    a_42
#         *      a_21   a_32    a_43 
#         a_11   a_22   a_33    A_44
#         a_12   a_23   a_44    *       ]
###


immutable BandedMatrix{T}
    data::Matrix{T}  # l+u+1 x n (# of rows)
    m::Int #Number of columns    
    l::Int # lower bandwidth ≥0
    u::Int # upper bandwidth ≥0
    function BandedMatrix(data::Matrix{T},m,l,u)
        @assert size(data,1)==l+u+1
        new(data,l,u,m)
    end    
end


BandedMatrix{T}(data::Matrix{T},m::Integer,a::Integer,b::Integer)=BandedMatrix{T}(data,m,a,b)
BandedMatrix{T}(data::Matrix{T},m::Integer,a)=BandedMatrix(data,m,-a[1],a[end])

BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a::Integer,b::Integer)=BandedMatrix{T}(Array(T,b+a+1,n),m,a,b)
BandedMatrix{T}(::Type{T},n::Integer,a::Integer,b::Integer)=BandedMatrix(T,n,n,a,b)
BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a)=BandedMatrix(T,n,m,-a[1],a[end])
BandedMatrix{T}(::Type{T},n::Integer,a)=BandedMatrix(T,n,-a[1],a[end])

Base.eltype{T}(::BandedMatrix{T})=T

for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::BandedMatrix,x::Number)=BandedMatrix($OP(B.data,x),B.m,B.l,B.u)
        $OP(x::Number,B::BandedMatrix)=BandedMatrix($OP(x,B.data),B.m,B.l,B.u)    
    end    
end




for (op,bop) in ((:(Base.rand),:barand),(:(Base.zeros),:bazeros),(:(Base.ones),:baones))
    @eval begin
        $bop{T}(::Type{T},n::Integer,m::Integer,a::Integer,b::Integer)=BandedMatrix($op(T,b+a+1,n),m,a,b)
        $bop{T}(::Type{T},n::Integer,m::Integer,a)=BandedMatrix($op(T,b+a+1,n),m,-a[1],a[end])        
        $bop{T}(::Type{T},n::Integer,a::Integer,b::Integer)=$bop(T,n,n,a,b)
        $bop{T}(::Type{T},n::Integer,a)=$bop(T,n,n,-a[1],a[end])        
        $bop(n::Integer,m::Integer,a::Integer,b::Integer)=$bop(Float64,n,m,a,b)
        $bop(n::Integer,m::Integer,a)=$bop(Float64,n,m,-a[1],a[end])        
        $bop(n::Integer,a::Integer,b::Integer)=$bop(n,n,a,b)
        $bop(n::Integer,a)=$bop(n,-a[1],a[end])        
    end
end



Base.size(A::BandedMatrix,k)=ifelse(k==1,size(A.data,2),A.m)
Base.size(A::BandedMatrix)=size(A.data,2),A.m
bandinds(A::BandedMatrix)=-A.l,A.u

usgetindex(A::BandedMatrix,k::Integer,j::Integer)=A.data[j-k+A.l+1,k]
usgetindex(A::BandedMatrix,k::Integer,jr::Range)=vec(A.data[jr-k+A.l+1,k])
getindex{T}(A::BandedMatrix{T},k::Integer,j::Integer)=(-A.l≤j-k≤A.u)?usgetindex(A,k,j):(j≤A.m?zero(T):throw(BoundsError()))
getindex(A::BandedMatrix,k::Integer,jr::Range)=-A.l≤jr[1]-k≤jr[end]-k≤A.u?usgetindex(A,k,jr):[A[k,j] for j=jr].'
getindex(A::BandedMatrix,k::Range,j::Integer)=[A[k,j] for j=jr]
getindex(A::BandedMatrix,kr::Range,jr::Range)=[A[k,j] for k=kr,j=jr]
Base.full(A::BandedMatrix)=A[1:size(A,1),1:size(A,2)]


# We turn off bound checking to allow nicer syntax without branching
#setindex!(A::BandedMatrix,v,k::Integer,j::Integer)=((A.l≤j-k≤A.u)&&k≤A.n)?ussetindex!(A,v,k,j):throw(BoundsError())
#setindex!(A::BandedMatrix,v,kr::Range,j::Integer)=(A.l≤j-kr[end]≤j-kr[1]≤A.u&&kr[end]≤A.n)?ussetindex!(A,v,kr,j):throw(BoundsError())


ibsetindex!(A::BandedMatrix,v,k::Integer,j::Integer)=(@inbounds A.data[j-k+A.l+1,k]=v)
ibpluseq!(A::BandedMatrix,v,k::Integer,j::Integer)=(@inbounds A.data[j-k+A.l+1,k]+=v)
setindex!(A::BandedMatrix,v,k::Integer,j::Integer)=(A.data[j-k+A.l+1,k]=v)

function setindex!(A::BandedMatrix,v,kr::Range,jr::Range)
    for j in jr
        A[kr,j]=slice(v,:,j)
    end
end
function setindex!(A::BandedMatrix,v,k::Integer,jr::Range)
    for j in jr
        A[k,j]=v[j]
    end
end



function bamultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix)   
    n,m=size(C)
    for k=1:n  # rows of C
        for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
            @inbounds Aj=A.data[l-k+A.l+1,k]
            
            shA=-l+B.l+1
            shB=-k+C.l+l-B.l
            @simd for j=max(1,k-C.l,l-B.l)+shA:min(B.u+l,n)+shA # columns of C/B
                @inbounds C.data[j+shB,k]+=Aj*B.data[j,l]
            end
        end
    end 
    C
end



function *{T,V}(A::BandedMatrix{T},B::BandedMatrix{V})
    if size(A,2)!=size(B,1)
        throw(DimensionMismatch("*"))
    end
    n=size(A,1)
    m=size(B,2)    
    bamultiply!(bazeros(promote_type(T,V),A.l+B.l,A.u+B.u,n,m),A,B)
end




## ShiftMatrix


immutable ShiftMatrix{T}
    data::Matrix{T}  # l+u+1 x n (# of rows)
    l::Int # lower bandwidth ≥0
    u::Int # upper bandwidth ≥0
    function ShiftMatrix(data::Matrix{T},l,u)
        @assert size(data,1)==l+u+1
        new(data,l,u)
    end    
end


ShiftMatrix{T}(data::Matrix{T},a::Integer,b::Integer)=ShiftMatrix{T}(data,a,b)
ShiftMatrix{T}(::Type{T},n::Integer,a::Integer,b::Integer)=ShiftMatrix{T}(Array(T,b+a+1,n),a,b)
ShiftMatrix{T}(::Type{T},n::Integer,a)=ShiftMatrix(T,n,-a[1],a[end])

Base.eltype{T}(::ShiftMatrix{T})=T

for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::ShiftMatrix,x::Number)=ShiftMatrix($OP(B.data,x),B.l,B.u)
        $OP(x::Number,B::ShiftMatrix)=ShiftMatrix($OP(x,B.data),B.l,B.u)    
    end    
end

for (op,bop) in ((:(Base.rand),:sarand),(:(Base.zeros),:sazeros),(:(Base.ones),:saones))
    @eval begin
        $bop{T}(::Type{T},n::Integer,a::Integer,b::Integer)=ShiftMatrix($op(T,b+a+1,n),a,b)
        $bop(n::Integer,a::Integer,b::Integer)=$bop(Float64,n,a,b)
        $bop(n::Integer,a)=$bop(Float64,n,-a[1],a[end])        
    end
end



Base.size(A::ShiftMatrix,k...)=size(A.data,k...)
columninds(A::ShiftMatrix)=-A.l,A.u
columnrange(A::ShiftMatrix)=-A.l:A.u

getindex(A::ShiftMatrix,k,j)=A.data[j+A.l+1,k]
Base.full(A::ShiftMatrix)=A[1:size(A,1),columnrange(A)]


# We turn off bound checking to allow nicer syntax without branching
#setindex!(A::ShiftMatrix,v,k::Integer,j::Integer)=((A.l≤j-k≤A.u)&&k≤A.n)?ussetindex!(A,v,k,j):throw(BoundsError())
#setindex!(A::ShiftMatrix,v,kr::Range,j::Integer)=(A.l≤j-kr[end]≤j-kr[1]≤A.u&&kr[end]≤A.n)?ussetindex!(A,v,kr,j):throw(BoundsError())


ibsetindex!(A::ShiftMatrix,v,k,j)=(@inbounds A.data[j+A.l+1,k]=v)
ibsetindex!(A::ShiftMatrix,v,k::Range,j::Range)=(@inbounds A.data[j+A.l+1,k]=v.')
ibpluseq!(A::ShiftMatrix,v,k,j)=(@inbounds A.data[j+A.l+1,k]+=v)
ibpluseq!(A::ShiftMatrix,v,k::Range,j::Range)=(@inbounds A.data[j+A.l+1,k]+=v.')
setindex!(A::ShiftMatrix,v,k,j)=(A.data[j+A.l+1,k]=v)
setindex!(A::ShiftMatrix,v,k::Range,j::Range)=(A.data[j+A.l+1,k]=v.')






## Convert

ShiftMatrix(B::BandedMatrix)=ShiftMatrix(B.data,B.l,B.u)
BandedMatrix(S::ShiftMatrix)=BandedMatrix(S.data,size(S,2)+S.u,S.l,S.u)
