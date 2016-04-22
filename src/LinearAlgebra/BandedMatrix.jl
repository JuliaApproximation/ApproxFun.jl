


function pad!(A::BandedMatrix,n,m)
    A.data=pad(A.data,size(A.data,1),n)
    A.m=m
    A
end

pad!(A::BandedMatrix,n,::Colon)=pad!(A,n,n+A.u)  # Default is to get all columns





## Used to scam addentries! into thinking we are somewhere else

immutable IndexStride{S,T} <: AbstractMatrix{T}
    matrix::S
    rowindex::Int
    colindex::Int
    rowstride::Int
    colstride::Int
end
function IndexStride{S<:BandedMatrix}(mat::S,ri::Int,ci::Int,rs::Int,cs::Int)
    # its no longer banded unless the strides match
    @assert rs==cs
    IndexStride{S,eltype(S)}(mat,ri,ci,rs,cs)
end
IndexStride(mat,ri,ci)=IndexStride(mat,ri,ci,1,1)


Base.size(S::IndexStride,k)=k==1?div(size(S.matrix,k)-S.rowindex,S.rowstride):div(size(S.matrix,k)-S.colindex,S.colstride)
Base.linearindexing{IS<:IndexStride}(::Type{IS})=Base.LinearSlow()

getindex(S::IndexStride,k,j)=S.matrix[S.rowstride*k+S.rowindex,S.colstride*j+S.colindex]
setindex!(S::IndexStride,x,k,j)=(S.matrix[S.rowstride*k+S.rowindex,S.colstride*j+S.colindex]=x)
unsafe_pluseq!(S::IndexStride,x,k,j)=unsafe_pluseq!(S.matrix,x,S.rowstride*k+S.rowindex,S.colstride*j+S.colindex)


isbeye(kr::Range)=IndexStride(beye(length(kr)),1-first(kr),1-first(kr))

function isbzeros{T}(::Type{T},kr::UnitRange,jr::UnitRange,l::Integer,u::Integer)
    shft=kr[1]-jr[1]
    IndexStride(bzeros(T,length(kr),length(jr),l-shft,u+shft),1-kr[1],1-jr[1])
end


# These view the operation as taking a slice of an infinite dimensional matrix
isbzeros{T}(::Type{T},kr::UnitRange,::Colon,l::Integer,u::Integer)=isbzeros(T,kr,max(1,kr[1]-l):kr[end]+u,l,u)
isbzeros{T}(::Type{T},kr::Colon,jr::UnitRange,l::Integer,u::Integer)=isbzeros(T,max(1,jr[1]-u):jr[end]+l,jr,l,u)

isbzeros{T}(::Type{T},rws::UnitRange,cols::UnitRange,bnds)=isbzeros(T,rws,cols,-bnds[1],bnds[end])
isbzeros{T}(::Type{T},kr::UnitRange,::Colon,bnds)=isbzeros(T,kr,:,-bnds[1],bnds[end])
isbzeros{T}(::Type{T},kr::Colon,jr::UnitRange,bnds)=isbzeros(T,:,jr,-bnds[1],bnds[end])

isbzeros(rw::Union{UnitRange,Colon},bnds...)=isbzeros(Float64,rw,bnds...)



for OP in (:*,:.*,:+,:.+,:-,:.-)
    @eval begin
        $OP{S}(B::IndexStride{S,Bool},x::Bool)=IndexStride($OP(B.matrix,x),B.rowindex,B.colindex,B.rowstride,B.colstride)
        $OP{S}(x::Bool,B::IndexStride{S,Bool})=IndexStride($OP(x,B.matrix),B.rowindex,B.colindex,B.rowstride,B.colstride)
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


## BandedMatrix operations


bandrange(S::IndexStride)=S.rowstride*bandrange(S.matrix)+S.rowindex-S.colindex
bandinds(S::IndexStride)=(S.rowstride*bandinds(S.matrix,1)+S.rowindex-S.colindex,S.rowstride*bandinds(S.matrix,2)+S.rowindex-S.colindex)
columnrange(A,row::Integer)=max(1,row+bandinds(A,1)):row+bandinds(A,2)


function eachbandedindex{ST<:BandedMatrix}(B::IndexStride{ST})
    map((k,j)->((k-S.rowindex)÷S.rowstride,(k-S.colindex)÷S.colstride),collect(eachbandedindex(B.matrix)))
end


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
unsafe_pluseq!(S::IndexSlice,x,k,j)=unsafe_pluseq!(S.matrix,x,div(k-S.rowindex,S.rowstride),div(j-S.colindex,S.colstride))


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
function setindex!(S::IndexTranspose,x,k,j)
    if S.firstrow≤k+j≤S.lastrow
        S.matrix[j,k]=x
    end
    x
end

#unsafe_pluseq!(S.matrix,x,k+j,-j)

function bmultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
   if rs==cs==1 && ri==ci==0 && size(C,1)==size(A,1) && size(C,2) == size(B,2)
        A_mul_B!(C,A,B)
    else
        error("Reimplemnt")
    end
end


function bmultiply!(C::IndexStride,A,B,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    bmultiply!(C.matrix,A,B,ri*C.rowstride+C.rowindex,ci*C.colstride+C.colindex,rs*C.rowstride,cs*C.colstride)
    C
end


function bmultiply!(C::IndexSlice,A,B,ri::Integer=0,ci::Integer=0,rs::Integer=1,cs::Integer=1)
    @assert rs==C.rowstride==cs==C.colstride
    @assert mod(ri-C.rowindex,rs)==mod(ci-C.colindex,cs)==0
    # div(rs*k+ri -C.rowindex,C.rowstride)
    # k+div(ri -C.rowindex,C.rowstride)
    bmultiply!(C.matrix,A,B,div(ri -C.rowindex,C.rowstride),div(ci -C.colindex,C.colstride))
    C
end




# Unoptimized but more readible version

# function bmultiply!(C::BandedMatrix,A::BandedMatrix,B::BandedMatrix)
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





## addentries!
# TODO: Remove

function addentries!(B::BandedMatrix,c::Number,A,kr::Range,::Colon)
    for (k,j) in eachbandedindex(B)
        A[k,j] += c*B.data[j-k+B.l+1,k]
    end

    A
end


function addentries!{ST<:BandedMatrix}(B::IndexStride{ST},c::Number,A,kr::Range,::Colon)
    for (k,j) in eachbandedindex(B)
        A[k,j]+=c*B[k,j]
    end

    A
end



addentries!(B::Union{BandedMatrix,IndexTranspose,IndexStride},A,kr::Range,::Colon)=addentries!(B,1,A,kr,:)
