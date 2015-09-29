

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
#         a_12   a_23   a_34    *       ]
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

BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a::Integer,b::Integer)=BandedMatrix{T}(Array(T,b+a+1,n),m,a,b)
BandedMatrix{T}(::Type{T},n::Integer,a::Integer,b::Integer)=BandedMatrix(T,n,n,a,b)
BandedMatrix{T}(::Type{T},n::Integer,::Colon,a::Integer,b::Integer)=BandedMatrix(T,n,n+b,a,b)


BandedMatrix{T}(data::Matrix{T},m::Integer,a)=BandedMatrix(data,m,-a[1],a[end])
BandedMatrix{T}(::Type{T},n::Integer,m::Integer,a)=BandedMatrix(T,n,m,-a[1],a[end])
BandedMatrix{T}(::Type{T},n::Integer,::Colon,a)=BandedMatrix(T,n,:,-a[1],a[end])
BandedMatrix{T}(::Type{T},n::Integer,a)=BandedMatrix(T,n,-a[1],a[end])

Base.eltype{T}(::BandedMatrix{T})=T
Base.convert{V}(::Type{BandedMatrix{V}},M::BandedMatrix) = BandedMatrix{V}(convert(Matrix{V},M.data),M.m,M.l,M.u)
function Base.convert{BM<:BandedMatrix}(::Type{BM},M::Matrix)
    ret=BandedMatrix(eltype(BM)==Any?eltype(M):promote_type(eltype(BM),eltype(M)),size(M,1),size(M,2),size(M,1)-1,size(M,2)-1)
    for k=1:size(M,1),j=1:size(M,2)
        ret[k,j]=M[k,j]
    end
    ret
end


Base.promote_rule{T,V}(::Type{BandedMatrix{T}},::Type{BandedMatrix{V}})=BandedMatrix{promote_type(T,V)}



for (op,bop) in ((:(Base.rand),:barand),(:(Base.zeros),:bazeros),(:(Base.ones),:baones))
    @eval begin
        $bop{T}(::Type{T},n::Integer,m::Integer,a::Integer,b::Integer)=BandedMatrix($op(T,b+a+1,n),m,a,b)
        $bop{T}(::Type{T},n::Integer,a::Integer,b::Integer)=$bop(T,n,n,a,b)
        $bop{T}(::Type{T},n::Integer,::Colon,a::Integer,b::Integer)=$bop(T,n,n+b,a,b)
        $bop{T}(::Type{T},::Colon,m::Integer,a::Integer,b::Integer)=$bop(T,m+a,m,a,b)
        $bop(n::Integer,m::Integer,a::Integer,b::Integer)=$bop(Float64,n,m,a,b)
        $bop(n::Integer,a::Integer,b::Integer)=$bop(n,n,a,b)

        $bop{T}(::Type{T},n::Integer,m::Integer,a)=$bop(T,n,m,-a[1],a[end])
        $bop{T}(::Type{T},n::Number,::Colon,a)=$bop(T,n,:,-a[1],a[end])
        $bop{T}(::Type{T},::Colon,m::Integer,a)=$bop(T,:,m,-a[1],a[end])
        $bop{T}(::Type{T},n::Integer,a)=$bop(T,n,-a[1],a[end])
        $bop(n::Integer,m::Integer,a)=$bop(Float64,n,m,-a[1],a[end])
        $bop(n::Integer,a)=$bop(n,-a[1],a[end])
    end
end



function baeye{T}(::Type{T},n::Integer,a...)
    ret=bazeros(T,n,a...)
    for k=1:n
         ret[k,k]=one(T)
    end
    ret
end
baeye{T}(::Type{T},n::Integer)=baeye(T,n,0,0)
baeye(n::Integer)=baeye(n,0,0)
baeye(n::Integer,a...)=baeye(Float64,n,a...)



Base.size(A::BandedMatrix,k)=ifelse(k==1,size(A.data,2),A.m)
Base.size(A::BandedMatrix)=size(A.data,2),A.m
bandinds(A::BandedMatrix)=-A.l,A.u
bandrange(A::BandedMatrix)=-A.l:A.u


usgetindex(A::BandedMatrix,k::Integer,j::Integer)=A.data[j-k+A.l+1,k]
usgetindex(A::BandedMatrix,k::Integer,jr::Range)=vec(A.data[jr-k+A.l+1,k])
getindex(A::BandedMatrix,k::Integer,j::Integer)=(-A.l≤j-k≤A.u)?usgetindex(A,k,j):(j≤A.m?zero(eltype(A)):throw(BoundsError()))
getindex(A::BandedMatrix,k::Integer,jr::Range)=-A.l≤jr[1]-k≤jr[end]-k≤A.u?usgetindex(A,k,jr):[A[k,j] for j=jr].'
getindex(A::BandedMatrix,kr::Range,j::Integer)=[A[k,j] for k=kr]
getindex(A::BandedMatrix,kr::Range,jr::Range)=[A[k,j] for k=kr,j=jr]
Base.full(A::BandedMatrix)=A[1:size(A,1),1:size(A,2)]


function Base.sparse(B::BandedMatrix)
    i=Array(Int,length(B.data));j=Array(Int,length(B.data))
    n,m=size(B.data)
    Bm=size(B,2)
    vb=copy(vec(B.data))
    for k=1:n,ℓ=1:m
        i[k+n*(ℓ-1)]=ℓ
        jj=k+ℓ-B.l-1
        if jj <1 || jj > Bm
            vb[k+n*(ℓ-1)] = 0
        end
        j[k+n*(ℓ-1)]=min(max(jj,1),Bm)
    end
    sparse(i,j,vb)
end




# pass standard routines to full matrix

Base.norm(B::BandedMatrix,opts...)=norm(full(B),opts...)


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

function *{T,V}(A::BandedMatrix{T},B::Matrix{V})
    if size(A,2)!=size(B,1)
        throw(DimensionMismatch("*"))
    end
    n,m=size(A,1),size(B,2)
    bamultiply!(zeros(promote_type(T,V),n,m),A,B)
end


function *{T,V}(A::BandedMatrix{T},b::Vector{V})
    if size(A,2)!=length(b)
        throw(DimensionMismatch("*"))
    end
    n=size(A,1)
    bamultiply!(zeros(promote_type(T,V),n),A,b)
end



function Base.transpose(B::BandedMatrix)
    Bt=bazeros(size(B,2),size(B,1),B.u,B.l)
    for k=1:size(B,1),j=max(1,k-B.l):min(size(B,2),k+B.u)
       Bt[j,k]=B[k,j]
    end
    Bt
end

function Base.ctranspose(B::BandedMatrix)
    Bt=bazeros(size(B,2),size(B,1),B.u,B.l)
    for k=1:size(B,1),j=max(1,k-B.l):min(size(B,2),k+B.u)
       Bt[j,k]=conj(B[k,j])
    end
    Bt
end



function pad!(A::BandedMatrix,n,m)
    A.data=pad(A.data,size(A.data,1),n)
    A.m=m
    A
end

pad!(A::BandedMatrix,n,::Colon)=pad!(A,n,n+A.u)  # Default is to get all columns




function Base.diag{T}(A::BandedMatrix{T})
    n=size(A,1)
    @assert n==size(A,2)

    vec(A.data[A.l+1,1:n])
end


## Used to scam addentries! into thinking we are somewhere else

immutable IndexStride{S}
    matrix::S
    rowindex::Int
    colindex::Int
    rowstride::Int
    colstride::Int
end
function IndexStride{S<:BandedMatrix}(mat::S,ri::Int,ci::Int,rs::Int,cs::Int)
    # its no longer banded unless the strides match
    @assert rs==cs
    IndexStride{S}(mat,ri,ci,rs,cs)
end
IndexStride(mat,ri,ci)=IndexStride(mat,ri,ci,1,1)

getindex(S::IndexStride,k,j)=S.matrix[S.rowstride*k+S.rowindex,S.colstride*j+S.colindex]
setindex!(S::IndexStride,x,k,j)=(S.matrix[S.rowstride*k+S.rowindex,S.colstride*j+S.colindex]=x)
ibpluseq!(S::IndexStride,x,k,j)=ibpluseq!(S.matrix,x,S.rowstride*k+S.rowindex,S.colstride*j+S.colindex)


# Odd rows/columns become real part and even become imag part
# S Should be a Real-valued matrix
immutable IndexReIm{S}
    matrix::S
end


getindex(S::IndexReIm,k,j)=S.matrix[2k-1,2j-1]+im*S.matrix[2k,2j]
function setindex!(S::IndexReIm,x,k,j)
    S.matrix[2k-1,2j-1]=real(x)
    S.matrix[2k,2j]=imag(x)
    x
end





isbaeye(kr::Range)=IndexStride(baeye(length(kr)),1-first(kr),1-first(kr))

function isbazeros{T}(::Type{T},kr::UnitRange,jr::UnitRange,l::Integer,u::Integer)
    shft=kr[1]-jr[1]
    IndexStride(bazeros(T,length(kr),length(jr),l-shft,u+shft),1-kr[1],1-jr[1])
end


# These view the operation as taking a slice of an infinite dimensional matrix
isbazeros{T}(::Type{T},kr::UnitRange,::Colon,l::Integer,u::Integer)=isbazeros(T,kr,max(1,kr[1]-l):kr[end]+u,l,u)
isbazeros{T}(::Type{T},kr::Colon,jr::UnitRange,l::Integer,u::Integer)=isbazeros(T,max(1,jr[1]-u):jr[end]+l,jr,l,u)

isbazeros{T}(::Type{T},rws::UnitRange,cols::UnitRange,bnds)=isbazeros(T,rws,cols,-bnds[1],bnds[end])
isbazeros{T}(::Type{T},kr::UnitRange,::Colon,bnds)=isbazeros(T,kr,:,-bnds[1],bnds[end])
isbazeros{T}(::Type{T},kr::Colon,jr::UnitRange,bnds)=isbazeros(T,:,jr,-bnds[1],bnds[end])

isbazeros(rw::Union{UnitRange,Colon},bnds...)=isbazeros(Float64,rw,bnds...)



for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::IndexStride,x::Number)=IndexStride($OP(B.matrix,x),B.rowindex,B.colindex,B.rowstride,B.colstride)
        $OP(x::Number,B::IndexStride)=IndexStride($OP(x,B.matrix),B.rowindex,B.colindex,B.rowstride,B.colstride)


        function $OP{ST<:BandedMatrix,SV<:BandedMatrix}(A::IndexStride{ST},B::IndexStride{SV})
            # TODO: General implementation
            @assert A.rowindex==B.rowindex==A.colindex==B.colindex
            @assert A.rowstride==B.rowstride==A.colstride==B.colstride==1

            AB=$OP(A.matrix,B.matrix)
            IndexStride(AB,A.rowindex,A.colindex)
        end
    end
end





bandrange(S::IndexStride)=S.rowstride*bandrange(S.matrix)+S.rowindex-S.colindex
bandinds(S::IndexStride)=(S.rowstride*bandinds(S.matrix,1)+S.rowindex-S.colindex,S.rowstride*bandinds(S.matrix,2)+S.rowindex-S.colindex)
columnrange(A,row::Integer)=max(1,row+bandinds(A,1)):row+bandinds(A,2)


## IndexSlice
#  divides by stride instead of multiply

immutable IndexSlice{S}
    matrix::S
    rowindex::Int
    colindex::Int
    rowstride::Int
    colstride::Int
end
function IndexSlice{S<:BandedMatrix}(mat::S,ri::Int,ci::Int,rs::Int,cs::Int)
    # its no longer banded unless the strides match
    @assert rs==cs
    @assert mod(ri-ci,rs)==0

    IndexSlice{S}(mat,ri,ci,rs,cs)
end
IndexSlice(mat,ri,ci)=IndexStride(mat,ri,ci,1,1)


#TODO: what if k and j are not in the stride?
#       this isn't an issue for now since we are using Slice
#       to make a small array look larger
getindex(S::IndexSlice,k,j)=S.matrix[div(k-S.rowindex,S.rowstride),div(j-S.colindex,S.colstride)]
setindex!(S::IndexSlice,x,k,j)=(S.matrix[div(k-S.rowindex,S.rowstride),div(j-S.colindex,S.colstride)]=x)
ibpluseq!(S::IndexSlice,x,k,j)=ibpluseq!(S.matrix,x,div(k-S.rowindex,S.rowstride),div(j-S.colindex,S.colstride))


for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP(B::IndexSlice,x::Number)=IndexStride($OP(B.matrix,x),B.rowindex,B.colindex,B.rowstride,B.colstride)
        $OP(x::Number,B::IndexSlice)=IndexStride($OP(x,B.matrix),B.rowindex,B.colindex,B.rowstride,B.colstride)
    end
end


# the following assume rowindex == colindex
bandinds(S::IndexSlice)=(div(bandinds(S.matrix,1)+S.colindex-S.rowindex,S.rowstride),div(bandinds(S.matrix,2)+S.colindex-S.rowindex,S.rowstride))

## Transpose indices

type IndexTranspose{S}
    matrix::S
    firstrow::Int   # These allow us to control which rows are toched
    lastrow::Int
end

IndexTranspose(mat,fl)=IndexTranspose(mat,fl[1],fl[end])


getindex(S::IndexTranspose,k,j)=S.matrix[j,k]
setindex!(S::IndexTranspose,x,k,j)=(S.matrix[j,k]=x)
function setindex!(S::IndexTranspose,x,k,j)
    if S.firstrow≤k+j≤S.lastrow
        S.matrix[j,k]=x
    end
    x
end

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




function bamultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    n=size(A,1);m=size(B,2)
    for k=1:n  # rows of C
        for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
            @inbounds Aj=A.data[l-k+A.l+1,k]


            #  A[k,j] == A.data[j-k+A.l+1,k]
            shB=-l+B.l+1
            ks=rs*k+ri
            shC=ci-ks+C.l+1
            @simd for j=max(1,l-B.l):min(B.u+l,m) # columns of C/B
                @inbounds C.data[cs*j+shC,ks]+=Aj*B.data[j+shB,l]
            end
        end
    end
    C
end

function bamultiply!(C::Matrix,A::BandedMatrix,B::Matrix,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    n=size(A,1);m=size(B,2)
    for k=1:n  # rows of C
        for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
            @inbounds Aj=A.data[l-k+A.l+1,k]

             @simd for j=1:m # columns of C/B
                 @inbounds C[rs*k+ri,cs*j+ci]+=Aj*B[l,j]
             end
        end
    end
    C
end


function bamultiply!(C::IndexStride,A,B,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    bamultiply!(C.matrix,A,B,ri*C.rowstride+C.rowindex,ci*C.colstride+C.colindex,rs*C.rowstride,cs*C.colstride)
    C
end


function bamultiply!(C::IndexSlice,A,B,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    @assert rs==C.rowstride==cs==C.colstride
    @assert mod(ri-C.rowindex,rs)==mod(ci-C.colindex,cs)==0
    # div(rs*k+ri -C.rowindex,C.rowstride)
    # k+div(ri -C.rowindex,C.rowstride)
    bamultiply!(C.matrix,A,B,div(ri -C.rowindex,C.rowstride),div(ci -C.colindex,C.colstride))
    C
end




# Unoptimized but more readible version

# function bamultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix)
#     n,m=size(C)
#     for k=1:n  # ROWS
#         for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
#             Aj=A[k,l]
#
#             for j=max(1,l-B.l):min(B.u+l,m) # columns of C/B
#                 C[rs*k+ri,cs*j+ci]+=Aj*B[l,j]
#             end
#         end
#     end
#     C
# end



## Matrix.*Matrix

function .*(A::BandedMatrix,B::BandedMatrix)
    @assert size(A,1)==size(B,1)&&size(A,2)==size(B,2)

    l=min(A.l,B.l);u=min(A.u,B.u)
    ret=BandedMatrix(promote_type(eltype(A),eltype(B)),size(A,1),size(A,2),l,u)

    for k=1:size(A,1),j=max(1,k-l):min(size(A,2),k+u)
        @inbounds ret[k,j]=A[k,j]*B[k,j]
    end
    ret
end


#implements fliplr(flipud(A))
function fliplrud(A::BandedMatrix)
    n,m=size(A)
    l=A.u+n-m
    u=A.l+m-n
    ret=BandedMatrix(eltype(A),n,m,l,u)
    for k=1:n,j=max(1,k-l):min(m,k+u)
        @inbounds ret[k,j]=A[n-k+1,m-j+1]
    end
    ret
end







## addentries!

function addentries!(B::BandedMatrix,c::Number,A,kr::Range,::Colon)
    for k=intersect(kr,1:size(B,1)),j=intersect(k+bandrange(B),1:size(B,2))
        A[k,j] += c*B.data[j-k+B.l+1,k]
    end

    A
end


function addentries!{ST<:BandedMatrix}(B::IndexStride{ST},c::Number,A,kr::Range,::Colon)
    for k=kr,j=k+bandrange(B)
        ibpluseq!(A,c*B[k,j],k,j)
    end

    A
end



addentries!(B::Union{BandedMatrix,IndexTranspose,IndexStride},A,kr::Range,::Colon)=addentries!(B,1,A,kr,:)
