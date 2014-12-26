

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


type BandedMatrix{T} <: AbstractSparseMatrix{T,Int}
    data::Matrix{T}  # l+u+1 x n (# of rows)
    m::Int #Number of columns    
    l::Int # lower bandwidth ≥0
    u::Int # upper bandwidth ≥0
    function BandedMatrix(data::Matrix{T},m,l,u)
        @assert size(data,1)==l+u+1
        new(data,m,l,u)
    end    
end


BandedMatrix{T}(data::Matrix{T},m::Integer,a::Integer,b::Integer)=BandedMatrix{T}(data,m,a,b)
BandedMatrix{T}(data::Matrix{T},m::Integer,a)=BandedMatrix(data,m,-a[1],a[end])

BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a::Integer,b::Integer)=BandedMatrix{T}(Array(T,b+a+1,n),m,a,b)
BandedMatrix{T}(::Type{T},n::Integer,a::Integer,b::Integer)=BandedMatrix(T,n,n,a,b)
BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a)=BandedMatrix(T,n,m,-a[1],a[end])
BandedMatrix{T}(::Type{T},n::Integer,a)=BandedMatrix(T,n,-a[1],a[end])

Base.eltype{T}(::BandedMatrix{T})=T





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
bandrange(A::BandedMatrix)=-A.l:A.u


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


for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::BandedMatrix{Bool},x::Bool)=BandedMatrix($OP(B.data,x),B.m,B.l,B.u)    
        $OP(x::Bool,B::BandedMatrix{Bool})=BandedMatrix($OP(x,B.data),B.m,B.l,B.u)            
        $OP(B::BandedMatrix,x::Number)=BandedMatrix($OP(B.data,x),B.m,B.l,B.u)
        $OP(x::Number,B::BandedMatrix)=BandedMatrix($OP(x,B.data),B.m,B.l,B.u)    
    end    
end

function +{T,V}(A::BandedMatrix{T},B::BandedMatrix{V})
    if size(A) != size(B)
        throw(DimensionMismatch("+"))
    end
    n,m=size(A,1),size(A,2)
    
    ret = bazeros(promote_type(T,V),n,m,max(A.l,B.l),max(A.u,B.u))
    for k=1:n,j=max(1,k-A.l):min(m,k+A.l)
        ibpluseq!(ret,usgetindex(A,k,j),k,j)
    end
    for k=1:n,j=max(1,k-B.l):min(m,k+B.l)
        ibpluseq!(ret,usgetindex(B,k,j),k,j)
    end
    
    ret    
end

function -{T,V}(A::BandedMatrix{T},B::BandedMatrix{V})
    if size(A) != size(B)
        throw(DimensionMismatch("+"))
    end
    n,m=size(A,1),size(A,2)
    
    ret = bazeros(promote_type(T,V),n,m,max(A.l,B.l),max(A.u,B.u))
    for k=1:n,j=max(1,k-A.l):min(m,k+A.l)
        ibpluseq!(ret,usgetindex(A,k,j),k,j)
    end
    for k=1:n,j=max(1,k-B.l):min(m,k+B.l)
        ibpluseq!(ret,-usgetindex(B,k,j),k,j)
    end
    
    ret    
end



function *{T,V}(A::BandedMatrix{T},B::BandedMatrix{V})
    if size(A,2)!=size(B,1)
        throw(DimensionMismatch("*"))
    end
    n,m=size(A,1),size(B,2)    
    bamultiply!(bazeros(promote_type(T,V),n,m,A.l+B.l,A.u+B.u),A,B)
end


function *{T,V}(A::BandedMatrix{T},b::Vector{V})
    if size(A,2)!=length(b)
        throw(DimensionMismatch("*"))
    end
    n=size(A,1)
    bamultiply!(zeros(promote_type(T,V),n),A,b)
end




## ShiftMatrix


type ShiftMatrix{T}
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
        $bop{T}(::Type{T},n::Integer,a)=$bop(T,n,-a[1],a[end])        
        $bop{T}(::Type{T},n::Integer)=$bop(T,n,0,0)          
        $bop(n::Integer,a::Integer,b::Integer)=$bop(Float64,n,a,b)
        $bop(n::Integer,a)=$bop(Float64,n,-a[1],a[end])        
        $bop(n::Integer)=$bop(n,0,0)          
    end
end

function saeye{T}(::Type{T},n::Integer,a...)
    ret=sazeros(T,n,a...)
    for k=1:n
         @inbounds ret[k,0]=one(T)
    end
    ret
end
saeye(n::Integer,a...)=saeye(Float64,n,a...)

baeye{T}(::Type{T},n::Integer,a...)=BandedMatrix(saeye(T,n,a...),n)
baeye(n::Integer,a...)=BandedMatrix(saeye(n,a...),n)


Base.size(A::ShiftMatrix,k)=ifelse(k==1,size(A.data,2),size(A.data,1))
Base.size(A::ShiftMatrix)=size(A,1),size(A,2)
columninds(A::ShiftMatrix)=-A.l,A.u
columnrange(A::ShiftMatrix)=-A.l:A.u
columninds(A::ShiftMatrix,k)=columninds(A)[k]

getindex(A::ShiftMatrix,k,j)=A.data[j+A.l+1,k]
Base.full(A::ShiftMatrix)=A[1:size(A,1),columnrange(A)]


# We turn off bound checking to allow nicer syntax without branching
#setindex!(A::ShiftMatrix,v,k::Integer,j::Integer)=((A.l≤j-k≤A.u)&&k≤A.n)?ussetindex!(A,v,k,j):throw(BoundsError())
#setindex!(A::ShiftMatrix,v,kr::Range,j::Integer)=(A.l≤j-kr[end]≤j-kr[1]≤A.u&&kr[end]≤A.n)?ussetindex!(A,v,kr,j):throw(BoundsError())


ibsetindex!(A::ShiftMatrix,v,k,j)=(@inbounds  A.data[j+A.l+1,k]=v)
ibsetindex!(A::ShiftMatrix,v,k::Range,j::Range)=(@inbounds  A.data[j+A.l+1,k]=v.')
ibpluseq!(A::ShiftMatrix,v,k,j)=(@inbounds  A.data[j+A.l+1,k]+=v)
ibpluseq!(A::ShiftMatrix,v,k::Range,j::Range)=(@inbounds  A.data[j+A.l+1,k]+=v.')
setindex!(A::ShiftMatrix,v,k,j)=(A.data[j+A.l+1,k]=v)
setindex!(A::ShiftMatrix,v,k::Range,j::Range)=(A.data[j+A.l+1,k]=v.')



function pad!(A::ShiftMatrix,n)
    A.data=pad(A.data,size(A.data,1),n)
    A
end



## Convert

ShiftMatrix(B::BandedMatrix)=ShiftMatrix(B.data,B.l,B.u)
BandedMatrix(S::ShiftMatrix)=BandedMatrix(S.data,size(S,1)+S.u,S.l,S.u)
BandedMatrix(S::ShiftMatrix,m)=BandedMatrix(S.data,m,S.l,S.u)



## Used to scam addentries! into thinking we are somewhere else

immutable IndexShift{S} 
    matrix::S
    rowindex::Int
    colindex::Int
    l::Int
    u::Int
end
IndexShift(S,ri,ci)=IndexShift(S,ri,ci,S.l+(ri-ci),S.u-(ri-ci))
IndexShift(S::ShiftMatrix,ri)=IndexShift(S,ri,0)
IndexShift(S::BandedMatrix,ri)=IndexShift(S,ri,ri)

getindex(S::IndexShift,k,j)=S.matrix[k-S.rowindex,j-S.colindex]
setindex!(S::IndexShift,x,k,j)=(S.matrix[k-S.rowindex,j-S.colindex]=x)
ibpluseq!(S::IndexShift,x,k,j)=ibpluseq!(S.matrix,x,k-S.rowindex,j-S.colindex)



columninds(S::IndexShift)=(columninds(S.matrix,1)+S.colindex,columninds(S.matrix,2)+S.colindex)
columnrange{BM<:BandedMatrix}(A::Union(IndexShift{BM},BandedMatrix),row::Integer)=max(1,row-A.l):row+A.u


for (isop,saop) in ((:issazeros,:sazeros),(:issaeye,:saeye))
    @eval begin
        $isop{T}(::Type{T},rws,bnds...)=IndexShift($saop(T,rws[end]-rws[1]+1,bnds...),rws[1]-1)
        $isop(rws,bnds...)=$isop(Float64,rws,bnds...)
    end
end

isbaeye(kr::Range)=IndexShift(baeye(length(kr)),first(kr)-1)


for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::IndexShift,x::Number)=IndexShift($OP(B.matrix,x),B.rowindex,B.colindex,B.l,B.u)
        $OP(x::Number,B::IndexShift)=IndexShift($OP(x,B.matrix),B.rowindex,B.colindex,B.l,B.u)

        
        function $OP{ST<:BandedMatrix,SV<:BandedMatrix}(A::IndexShift{ST},B::IndexShift{SV})
            # TODO: General implementation
            @assert A.rowindex==B.rowindex==A.colindex==B.colindex
            
            AB=$OP(A.matrix,B.matrix)
            IndexShift(AB,A.rowindex,A.colindex)
        end        
    end    
end





columnrange(S::IndexShift)=columnrange(S.matrix)+S.colindex
columninds(S::IndexShift)=(columninds(S.matrix,1)+S.colindex,columninds(S.matrix,2)+S.colindex)
bandrange(S::IndexShift)=bandrange(S.matrix)+S.colindex-S.rowindex
bandinds(S::IndexShift)=(bandinds(S.matrix,1)+S.colindex-S.rowindex,bandinds(S.matrix,2)+S.colindex-S.rowindex)


type IndexTranspose{S}
    matrix::S
    firstrow::Int   # These allow us to control which rows are toched
    lastrow::Int    
end

IndexTranspose(mat,fl)=IndexTranspose(mat,fl[1],fl[end])

getindex{ST<:ShiftMatrix}(S::IndexTranspose{ST},k,j)=S.matrix[k+j,-j]
function setindex!{ST<:ShiftMatrix}(S::IndexTranspose{ST},x,k,j)
    if S.firstrow≤k+j≤S.lastrow
        S.matrix[k+j,-j]=x
    end
    x
end

getindex(S::IndexTranspose,k,j)=S.matrix[j,k]
setindex!(S::IndexTranspose,x,k,j)=(S.matrix[j,k]=x)
function setindex!(S::IndexTranspose,x,k,j)
    if S.firstrow≤k+j≤S.lastrow
        S.matrix[j,k]=x
    end
    x
end
ibpluseq!{ST<:ShiftMatrix}(S::IndexTranspose{ST},x,k,j)=S.matrix[k+j,-j]+=x
#ibpluseq!(S.matrix,x,k+j,-j)

## Matrix*Vector Multiplicaiton

function bamultiply!(c::Vector,A::BandedMatrix,b::Vector)
    for k=1:size(A,1)  # rows of c
        @simd for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A/rows of b
             @inbounds c[k]+=A.data[l-k+A.l+1,k]*b[l]
        end        
    end
    c
end




## Matrix*Matrix Multiplication

function bamultiply!(C::Union(BandedMatrix,ShiftMatrix),A::BandedMatrix,B::BandedMatrix,rs::Integer=0,cs::Integer=0)   
    n=size(A,1);m=size(B,2)
    for k=1:n  # rows of C
        for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
            @inbounds Aj=A.data[l-k+A.l+1,k]
            
            shA=-l+B.l+1
            shB=-k+C.l+l-B.l+cs-rs
            @simd for j=(max(1-shB,k-C.l,l-B.l)+shA):(min(B.u+l,m)+shA) # columns of C/B
                @inbounds C.data[j+shB,k+rs]+=Aj*B.data[j,l]
            end
        end
    end 
    C
end

function bamultiply!{SM<:ShiftMatrix}(C::IndexShift{SM},A::BandedMatrix,B::BandedMatrix,rs::Integer=0,cs::Integer=0)
    bamultiply!(C.matrix,A,B,rs-C.rowindex,cs-C.colindex-C.rowindex)
    C
end


function bamultiply!{SM<:BandedMatrix}(C::IndexShift{SM},A::BandedMatrix,B::BandedMatrix,rs::Integer=0,cs::Integer=0)
    bamultiply!(C.matrix,A,B,rs-C.rowindex,cs-C.colindex)
    C
end




# Unoptimized but more readible version

# function bamultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix)   
#     n,m=size(C)
#     for k=1:n  # ROWS
#         for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
#             Aj=A[k,l]
# 
#             for j=max(1,k-C.l,l-B.l):min(B.u+l,n) # columns of C/B
#                 C[k,j]+=Aj*B[l,j]
#             end
#         end
#     end 
#     C
# end






## addentries!

function addentries!(B::BandedMatrix,c::Number,A,kr::Range)    
    for k=kr,j=k+bandrange(B)
        A[k,j-k] += c*B[k,j]
    end
    
    A
end

function addentries!(B::ShiftMatrix,c::Number,A,kr::Range)    
    for k=kr,j=columnrange(B)
         @inbounds ibpluseq!(A,c*B.data[j+B.l+1,k],k,j)
    end
    
    A
end

function addentries!{ST<:ShiftMatrix}(B::IndexShift{ST},c::Number,A,kr::Range)    
    for k=kr,j=max(1-k,columninds(B,1)):columninds(B,2)
        ibpluseq!(A,c*B[k,j],k,j)
    end
    
    A
end

function addentries!{ST<:BandedMatrix}(B::IndexShift{ST},c::Number,A,kr::Range)    
    for k=kr,j=bandrange(B)
        ibpluseq!(A,c*B[k,k+j],k,j)
    end
    
    A
end


# function addentries!{ST<:ShiftMatrix}(B::IndexTranpose{ST},c::Number,A,kr::Range)    
#     for k=kr,j=columnrange(B)
#         A[k,j] += c*B[k,j]
#     end
#     
#     A
# end



addentries!(B::Union(BandedMatrix,ShiftMatrix,IndexTranspose,IndexShift),A,kr::Range)=addentries!(B,1,A,kr)
