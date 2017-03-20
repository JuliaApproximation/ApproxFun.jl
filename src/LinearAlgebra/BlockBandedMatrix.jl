doc"""
`Block` is used for get index of a block matrix.  For example,
```julia
A[Block(1),Block(2)]
```
retrieves the 1 x 2 block
"""
immutable Block <: Real
    K::Int
end

block(B::Block) = B

Base.convert(::Type{Integer},B::Block) = B.K
Base.convert{T<:Integer}(::Type{T},B::Block) = convert(T,B.K)::T
Base.convert(::Type{Block},K::Block) = K
Base.convert(::Type{Block},K::Integer) = Block(K)

Base.promote_rule(::Type{Block},::Type{Int}) = Block

# supports SubArray
Base.to_index(b::Block) = b

for OP in (:(Base.one),:(Base.zero),:(-),:(Base.floor))
    @eval $OP(B::Block) = Block($OP(B.K))
end

for OP in (:(==),:(!=),:(Base.isless),:(<=),:(>=),:(<),:(>))
    @eval $OP(A::Block,B::Block) = $OP(A.K,B.K)
end

for OP in (:(Base.iseven),:(Base.isodd))
    @eval $OP(A::Block) = $OP(A.K)
end

Base.isless(::Block,::Infinity{Bool}) = true
Base.isless(::Infinity{Bool},::Block) = false


for OP in (:(Base.rem),:(Base.div))
    @eval begin
        $OP(A::Block,B::Block) = Block($OP(A.K,B.K))
        $OP(A::Integer,B::Block) = Block($OP(A,B.K))
        $OP(A::Block,B::Integer) = Block($OP(A.K,B))
    end
end


# broadcasting
Base.broadcast(::Type{Int},r::UnitRange{Block}) = Int(r.start):Int(r.stop)
Base.broadcast(::Type{Int},r::StepRange{Block,Block}) = Int(r.start):Int(r.step):Int(r.stop)
Base.broadcast(::Type{Int},B::Block) = Int(B)
#
Base.broadcast(::Type{Block},r::UnitRange{Int}) = Block(r.start):Block(r.stop)
Base.broadcast(::Type{Block},r::StepRange{Int,Int}) = Block(r.start):Block(r.step):Block(r.stop)
Base.broadcast(::Type{Block},B::Int) = Block(B)



doc"""
`SubBlock` is used for get subindices of a block.  For example,
```julia
A[Block(1)[1:3],Block(2)[3:4]]
```
retrieves 1:3 × 3:4 entries of the 1 × 2 block
"""
immutable SubBlock{R}
    block::Block
    sub::R
end

block(B::SubBlock) = B.block
getindex(A::Block,k) = SubBlock(A,k)
getindex(A::SubBlock,k) = SubBlock(A.block,A.sub[k])



# supports SubArray
Base.to_index(b::SubBlock) = b
Base.size{RR<:Range}(b::SubBlock{RR}) = size(b.sub)


for OP in (:blockrows, :blockcols, :blockrowstart, :blockrowstop, :blockcolstart, :blockcolstop)
    @eval $OP(A,K::Block) = $OP(A,K.K)
end

for OP in (:+,:-,:*)
    @eval begin
        $OP(K::Block,J::Integer) = Block($OP(K.K,J))
        $OP(J::Integer,K::Block) = Block($OP(J,K.K))
        $OP(K::Block,J::Block) = Block($OP(K.K,J.K))
    end
end

# Takes in a list of lengths, and allows lookup of which block an entry is in
function blocklookup(rows)
    rowblocks=Array{Int}(0)
    for ν in eachindex(rows), k in 1:rows[ν]
        push!(rowblocks,ν)
    end
    rowblocks
end

if VERSION ≥ v"0.6-"
    # 0.6 index routines
    @inline Base.index_dimsum(::Block, I...) = (true, Base.index_dimsum(I...)...)
    @inline Base.index_dimsum(::SubBlock, I...) = (true, Base.index_dimsum(I...)...)
end


# AbstractBlockMatrix`
@compat abstract type AbstractBlockMatrix{T,BT} <: AbstractMatrix{T} end


rowblocklengths(A::AbstractBlockMatrix) = A.rows
colblocklengths(A::AbstractBlockMatrix) = A.cols


function getindex{T,BT}(A::AbstractBlockMatrix{T,BT},K::Block,J::Block)
    if -A.l ≤ J-K ≤ A.u
        BT(view(A,K,J))
    else
        zeroblock(A,K,J)
    end
end
getindex(A::AbstractBlockMatrix,K::Block,j) = A[blockrows(A,K),j]
getindex(A::AbstractBlockMatrix,k,J::Block) = A[k,blockcols(A,J)]

getindex(A::AbstractBlockMatrix,KR::UnitRange{Block},JR::UnitRange{Block}) =
    Matrix(view(A,KR,JR))

setindex!(A::AbstractBlockMatrix,V,K::Block,J::Block) = (view(A,K,J) .= V)
setindex!(A::AbstractBlockMatrix,V,K::Block,j) = (view(A,blockrows(A,K),j) .= V)
setindex!(A::AbstractBlockMatrix,V,k,J::Block) = (view(A,k,blockcols(A,J)) .= V)


Base.size(A::AbstractBlockMatrix) = sum(A.rows),sum(A.cols)

blocksize(A::AbstractMatrix,k::Int) = k==1?length(rowblocklengths(A)):length(colblocklengths(A))
blocksize(A::AbstractMatrix) = blocksize(A,1),blocksize(A,2)


function blockmatrix_αA_mul_B_plus_βC!(α,A,x,β,y)
    if length(x) != size(A,2) || length(y) != size(A,1)
        throw(BoundsError())
    end

    BLAS.scal!(length(y),β,y,1)
    o=one(eltype(y))

    for J=Block(1):Block(blocksize(A,2))
        jr=blockcols(A,J)
        for K=blockcolrange(A,J)
            kr=blockrows(A,K)
            B=view(A,K,J)
            αA_mul_B_plus_βC!(α,B,view(x,jr),o,view(y,kr))
        end
    end
    y
end

αA_mul_B_plus_βC!(α,A::AbstractBlockMatrix,x::AbstractVector,β,y::AbstractVector) =
    blockmatrix_αA_mul_B_plus_βC!(α,A,x,β,y)

Base.A_mul_B!(y::Vector,A::AbstractBlockMatrix,b::Vector) =
    αA_mul_B_plus_βC!(one(eltype(A)),A,b,zero(eltype(y)),y)




function block_axpy!(α,A,Y)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=Block.(1:blocksize(A,2)), K=blockcolrange(A,J)
        BLAS.axpy!(α,view(A,K,J),view(Y,K,J))
    end
    Y
end

Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::AbstractBlockMatrix) = block_axpy!(α,A,Y)



## Algebra


function Base.A_mul_B!(Y::AbstractBlockMatrix,A::AbstractBlockMatrix,B::AbstractBlockMatrix)
    T=eltype(Y)
    BLAS.scal!(length(Y.data),zero(T),Y.data,1)
    o=one(T)
    for J=Block(1):Block(blocksize(B,2)),N=blockcolrange(B,J),K=blockcolrange(A,N)
        αA_mul_B_plus_βC!(o,view(A,K,N),view(B,N,J),o,view(Y,K,J))
    end
    Y
end

## Sub Block

blocksize{T,BBM<:AbstractBlockMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}},k::Int) =
    k == 1 ? parent(S).rowblocks[last(parentindexes(S)[1])] : parent(S).colblocks[last(parentindexes(S)[2])]

function rowblocklengths{T,BBM<:AbstractBlockMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}})
    P = parent(S)
    kr = parentindexes(S)[1]
    B1=P.rowblocks[kr[1]]
    B2=P.rowblocks[kr[end]]

    ret = zeros(Int,B2)
    # if the blocks are equal, we have only one bvlock
    if B1 == B2
        ret[end] = length(kr)
    else
        ret[B1] = blockrows(P,B1)[end]-kr[1]+1
        ret[B1+1:B2-1] = view(P.rows,B1+1:B2-1)
        ret[end] = kr[end]-blockrows(P,B2)[1]+1
    end
    ret
end

function colblocklengths{T,BBM<:AbstractBlockMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}})
    P = parent(S)
    jr = parentindexes(S)[2]
    B1=P.colblocks[jr[1]]
    B2=P.colblocks[jr[end]]

    ret = zeros(Int,B2)
    # if the blocks are equal, we have only one bvlock
    if B1 == B2
        ret[end] = length(jr)
    else
        ret[B1] = blockcols(P,B1)[end]-jr[1]+1
        ret[B1+1:B2-1] = view(P.cols,B1+1:B2-1)
        ret[end] = jr[end]-blockcols(P,B2)[1]+1
    end
    ret
end


function blockrows{T,BBM<:AbstractBlockMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}},K::Int)
    kr = parentindexes(S)[1]
    br = blockrows(parent(S),K)-kr[1]+1
    max(1,br[1]):min(length(kr),br[end])
end

function blockcols{T,BBM<:AbstractBlockMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}},J::Int)
    jr = parentindexes(S)[2]
    br = blockcols(parent(S),J)-jr[1]+1
    max(1,br[1]):min(length(jr),br[end])
end

@compat const AllBlockMatrix{T,BBM<:AbstractBlockMatrix,II<:Union{UnitRange{Int},UnitRange{Block}},JJ<:Union{UnitRange{Int},UnitRange{Block}}} =
    Union{AbstractBlockMatrix,SubArray{T,2,BBM,Tuple{II,JJ}}}

function Base.BLAS.axpy!(α,A::AllBlockMatrix,Y::AbstractMatrix)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=1:blocksize(A,2)
        jr=blockcols(A,J)
        for K=blockcolrange(A,J)
            kr=blockrows(A,K)
            BLAS.axpy!(α,view(A,Block(K),Block(J)),view(Y,kr,jr))
        end
    end
    Y
end


Base.BLAS.axpy!(α,A::AllBlockMatrix,Y::AllBlockMatrix) = block_axpy!(α,A,Y)




## AbstractBlockBandedMatrix

@compat abstract type AbstractBlockBandedMatrix{T,BT} <: AbstractBlockMatrix{T,BT} end

isblockbanded(::) = false
isblockbanded(::AbstractBlockBandedMatrix) = true


blockbandinds(K::AbstractBlockBandedMatrix) = (-K.l,K.u)
blockbandwidths(S) = (-blockbandinds(S,1),blockbandinds(S,2))
blockbandinds(K,k::Integer) = blockbandinds(K)[k]
blockbandwidth(K,k::Integer) = k==1 ? -blockbandinds(K,k) : blockbandinds(K,k)

# these give the block rows corresponding to the J-th column block
blockcolstart(A::AbstractBlockBandedMatrix,J::Int) = Block(max(1,J-A.u))
blockcolstop(A::AbstractBlockBandedMatrix,J::Int) = Block(min(length(A.rows),J+A.l))
blockcolstart(A::AbstractMatrix,J::Int) = Block(max(1,J-blockbandwidth(A,2)))
blockcolstop(A::AbstractMatrix,J::Int) = Block(min(blocksize(A,1),J+blockbandwidth(A,1)))
blockcolrange(A,J) = blockcolstart(A,J):blockcolstop(A,J)

blockrowstart(A::AbstractBlockBandedMatrix,K::Int) = Block(max(1,K-A.l))
blockrowstop(A::AbstractBlockBandedMatrix,K::Int) = Block(min(length(A.cols),K+A.u))
blockrowstart(A::AbstractMatrix,K::Int) = Block(max(1,K-blockbandwidth(A,1)))
blockrowstop(A::AbstractMatrix,K::Int) = Block(min(length(colblocklengths(A)),K+blockbandwidth(A,2)))
blockrowrange(A,K) = blockrowstart(A,K):blockrowstop(A,K)

# give the rows/columns of a block, as a range
blockrows(A::AbstractBlockMatrix,K::Int) = sum(A.rows[1:K-1]) + (1:A.rows[K])
blockcols(A::AbstractBlockMatrix,J::Int) = sum(A.cols[1:J-1]) + (1:A.cols[J])

for Op in (:blockcolstart,:blockcolstop,:blockrowstart,:blockrowstop,:blockrows,:blockcols)
    @eval $Op(A::AbstractBlockMatrix,K::Block) = $Op(A,K.K)
end

colstop(A::AbstractBlockBandedMatrix,k::Int) = sum(A.rows[1:min(A.colblocks[k]+A.l,length(A.rows))])
function colstart(A::AbstractBlockBandedMatrix,k::Int)
    lr = A.colblocks[k]-A.u-1
    lr ≤ length(A.rows) ? sum(A.rows[1:lr])+1 : size(A,1)+1  # no columns
end

rowstop(A::AbstractBlockBandedMatrix,k::Int) = sum(A.cols[1:min(A.rowblocks[k]+A.u,length(A.cols))])
rowstart(A::AbstractBlockBandedMatrix,k::Int) = sum(A.cols[1:A.rowblocks[k]-A.l-1])+1

Base.convert(::Type{Matrix},A::AbstractBlockMatrix) =
    BLAS.axpy!(one(eltype(A)),A,zeros(eltype(A),size(A,1),size(A,2)))

Base.full(S::AbstractBlockMatrix) = convert(Matrix, S)


function Base.getindex(A::AbstractBlockBandedMatrix,k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        zero(eltype(A))
    else
        k2=k-sum(A.rows[1:K-1])
        j2=j-sum(A.cols[1:J-1])
        view(A,Block(K),Block(J))[k2,j2]
    end
end

function Base.setindex!(A::AbstractBlockBandedMatrix,v,k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        if v == 0
            return v
        else
            error("block $K, $J outside bands.")
        end
    else
        k2=k-sum(A.rows[1:K-1])
        j2=j-sum(A.cols[1:J-1])
        view(A,Block(K),Block(J))[k2,j2] = v
    end
end

Base.getindex(A::AbstractBlockBandedMatrix,kj::CartesianIndex{2}) =
    A[kj[1],kj[2]]

@compat Base.IndexStyle{BBBM<:AbstractBlockMatrix}(::Type{BBBM}) =
    IndexCartesian()




## BlockBandedMatrix routines


function bbm_blockstarts(rows,cols,l,u)
    numentries = 0
    N=length(rows)
    # TODO: remove max when BandedMatrix allows negative indices
    bs=BandedMatrix(Int,length(rows),length(cols),max(0,l),max(0,u))
    for J = eachindex(cols), K = max(1,J-u):min(J+l,N)
        bs[K,J] = numentries
        numentries += cols[J]*rows[K]
    end
    bs
end


function bbm_numentries(rows,cols,l,u)
    numentries = 0
    N=length(rows)
    for J = eachindex(cols), K = max(1,J-u):min(J+l,N)
        numentries += cols[J]*rows[K]
    end
    numentries
end


#  A block matrix where only theb ands are nonzero
#   isomorphic to BandedMatrix{Matrix{T}}
type BlockBandedMatrix{T,RI,CI} <: AbstractBlockBandedMatrix{T,Matrix{T}}
    data::Vector{T}

    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    rows::RI   # the size of the row blocks
    cols::CI  # the size of the column blocks

    rowblocks::Vector{Int}   # lookup the row block given the row
    colblocks::Vector{Int}   # lookup col row block given the row

    blockstart::BandedMatrix{Int}  # gives shift of first entry of data for each block

    function (::Type{BlockBandedMatrix{T,RI,CI}}){T,RI,CI}(data::Vector{T},l,u,rows,cols)
        new{T,RI,CI}(data,l,u,rows,cols,blocklookup(rows),blocklookup(cols),bbm_blockstarts(rows,cols,l,u))
    end
end






BlockBandedMatrix(data::Vector,l,u,rows,cols) =
    BlockBandedMatrix{eltype(data),typeof(rows),typeof(cols)}(data,l,u,rows,cols)

BlockBandedMatrix{T}(::Type{T},l,u,rows,cols) =
    BlockBandedMatrix(Array{T}(bbm_numentries(rows,cols,l,u)),l,u,rows,cols)

for FUNC in (:zeros,:rand,:ones)
    BFUNC = parse("bb"*string(FUNC))
    @eval $BFUNC{T}(::Type{T},l,u,rows,cols) =
        BlockBandedMatrix($FUNC(T,bbm_numentries(rows,cols,l,u)),l,u,rows,cols)
end

#TODO: override and use Base.copy!
function Base.convert(::Type{BlockBandedMatrix},Y::AbstractBlockMatrix)
    ret = bbzeros(eltype(Y),Y.l,Y.u,Y.rows,Y.cols)
    BLAS.axpy!(one(eltype(Y)),Y,ret)
    ret
end

function Base.convert(::Type{BlockBandedMatrix},Y::AbstractMatrix)
    if !isblockbanded(Y)
        error("Cannot convert $(typeof(Y)) is not block banded")
    end
    ret = bbzeros(eltype(Y),blockbandwidth(Y,1),blockbandwidth(Y,2),rowblocklengths(Y),colblocklengths(Y))
    BLAS.axpy!(one(eltype(Y)),Y,ret)
    ret
end


zeroblock(X::BlockBandedMatrix,K::Block,J::Block) = zeros(eltype(X),X.rows[K.K],X.cols[J.K])



## View

@compat const SubBandedBlockSubBlock{T,BBM<:BlockBandedMatrix,II<:Union{Block,SubBlock},JJ<:Union{Block,SubBlock}} = SubArray{T,2,BBM,Tuple{II,JJ},false}
@compat const SubBandedBlockRange{T,BBM<:BlockBandedMatrix} = SubArray{T,2,BBM,Tuple{UnitRange{Block},UnitRange{Block}},false}
@compat const StridedMatrix2{T,A<:Union{DenseArray,Base.StridedReshapedArray,BlockBandedMatrix},I<:Tuple{Vararg{Union{Base.RangeIndex,Base.AbstractCartesianIndex,Block,SubBlock}}}} = Union{DenseArray{T,2}, SubArray{T,2,A,I}, Base.StridedReshapedArray{T,2}}


isblockbanded(::SubBandedBlockRange) = true
isblockbanded{T,BBM<:AbstractBlockBandedMatrix}(::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}},false}) = true

subblocksize(A::AbstractBlockBandedMatrix,K::Block,J::Block)::Tuple{Int,Int} =
    (A.rows[K.K],A.cols[J.K])
subblocksize(A::AbstractBlockBandedMatrix,K::Block,J::SubBlock)::Tuple{Int,Int} =
    (A.rows[K.K],length(J.sub))
subblocksize(A::AbstractBlockBandedMatrix,K::SubBlock,J::Block)::Tuple{Int,Int} =
    (length(K.sub),A.cols[J.K])
subblocksize(A::AbstractBlockBandedMatrix,K::SubBlock,J::SubBlock)::Tuple{Int,Int} =
    (length(K.sub),length(J.sub))

subblocksize(A::AbstractBlockBandedMatrix,K::UnitRange{Block},J::UnitRange{Block})::Tuple{Int,Int} =
    (sum(A.rows[Int.(K)]),sum(A.cols[Int.(J)]))

rowblocklengths(A::SubBandedBlockRange) = parent(A).rows[parentindexes(A)[1]]
colblocklengths(A::SubBandedBlockRange) = parent(A).cols[parentindexes(A)[2]]


 ## Indices
function Base.indices(S::Union{SubBandedBlockSubBlock,SubBandedBlockRange})
    sz = subblocksize(parent(S),parentindexes(S)...)::Tuple{Int,Int}
    (Base.OneTo(sz[1]),
     Base.OneTo(sz[2]))
 end
 Base.indices{T,BB<:AbstractBlockBandedMatrix}(S::SubArray{T,2,BB,Tuple{Block,Block},false}) =
     (Base.OneTo(Int(parent(S).rows[parentindexes(S)[1].K])),
      Base.OneTo(Int(parent(S).cols[parentindexes(S)[2].K])))
Base.indices{T,BB<:AbstractBlockBandedMatrix}(S::SubArray{T,2,BB,Tuple{SubBlock{UnitRange{Int}},SubBlock{UnitRange{Int}}},false}) =
  (Base.OneTo(length(parentindexes(S)[1].sub)::Int),
   Base.OneTo(length(parentindexes(S)[2].sub)::Int))


parentblock(S) = view(parent(S),block(parentindexes(S)[1]),block(parentindexes(S)[2]))


# returns a view of the data corresponding to the block
function blockview{T,U,V}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{Block,Block},false})
    A = parent(S)
    K,J = parentindexes(S)
    st = A.blockstart[block(K).K,block(J).K]
    sz = size(S)
    reshape(view(A.data,st+1:st+sz[1]*sz[2]),sz[1],sz[2])
end

# returns a view of the data
dataview{T,U,V}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{Block,Block},false}) =
    blockview(S)
dataview{T,U,V,II}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{SubBlock{II},Block},false}) =
    view(blockview(parentblock(S)),parentindexes(S)[1].sub,:)
dataview{T,U,V,JJ}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{Block,SubBlock{JJ}},false}) =
    view(blockview(parentblock(S)),:,parentindexes(S)[2].sub)
dataview{T,U,V,II,JJ}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{SubBlock{II},SubBlock{JJ}},false}) =
    view(blockview(parentblock(S)),parentindexes(S)[1].sub,parentindexes(S)[2].sub)


blockbandinds{T,BBM<:AbstractBlockBandedMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}}) =
    blockbandinds(parent(S))

getindex{N}(S::SubBandedBlockSubBlock, I::Vararg{Int,N}) =
    dataview(S)[I...]

function getindex(S::SubBandedBlockRange, k::Integer, j::Integer)
    KR,JR = parentindexes(S)
    A = parent(S)
    A[k + sum(A.rows[1:first(KR).K-1]),j + sum(A.rows[1:first(JR).K-1])]
end

function setindex!{N}(S::SubBandedBlockSubBlock, v, I::Vararg{Int,N})
    dataview(S)[I...] = v
end

function setindex!(S::SubBandedBlockRange, v, k::Integer, j::Integer)
    KR,JR = parentindexes(S)
    A = parent(S)
    A[k + sum(A.rows[1:first(KR).K-1]),j + sum(A.rows[1:first(JR).K-1])] = v
end

Base.BLAS.axpy!(α,A::SubBandedBlockRange,Y::SubBandedBlockRange) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::SubBandedBlockRange,Y::AbstractBlockMatrix) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::SubBandedBlockRange) = block_axpy!(α,A,Y)


Base.strides(S::SubBandedBlockSubBlock) =
    (1,parent(S).rows[block(parentindexes(S)[1]).K])

function Base.pointer{T<:BlasFloat,U,V}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{Block,Block},false})
    A = parent(S)
    K,J = parentindexes(S)
    pointer(A.data)+A.blockstart[K.K,J.K]*sizeof(T)
end


Base.pointer{T<:BlasFloat,U,V,JJ}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{Block,SubBlock{JJ}},false}) =
    pointer(parentblock(S))

Base.pointer{T<:BlasFloat,U,V,II}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{SubBlock{II},Block},false}) =
    pointer(parentblock(S)) + sizeof(T)*(first(parentindexes(S)[1].sub)-1)

Base.pointer{T<:BlasFloat,U,V,II,JJ}(S::SubArray{T,2,BlockBandedMatrix{T,U,V},Tuple{SubBlock{II},SubBlock{JJ}},false}) =
    pointer(parentblock(S)) + sizeof(T)*(first(parentindexes(S)[1].sub)-1) +
        sizeof(T)*stride(S,2)*(first(parentindexes(S)[2].sub)-1)


## algebra

αA_mul_B_plus_βC!{T}(α::T,A::SubBandedBlockSubBlock{T},x::AbstractVector{T},β::T,y::AbstractVector{T}) =
    gemv!('N',α,A,x,β,y)
αA_mul_B_plus_βC!{T}(α,A::SubBandedBlockSubBlock{T},x::AbstractVector,β,y::AbstractVector{T}) =
    gemv!('N',T(α),A,Array{T}(x),T(β),y)
αA_mul_B_plus_βC!(α,A::StridedMatrix2,B::StridedMatrix2,β,C::StridedMatrix2) =
    gemm!('N','N',α,A,B,β,C)

function *{T<:Number,V<:Number}(A::BlockBandedMatrix{T},
                                B::BlockBandedMatrix{V})
    if A.cols ≠ B.rows
        throw(DimensionMismatch("*"))
    end
    n,m=size(A,1),size(B,2)

    A_mul_B!(BlockBandedMatrix(promote_type(T,V),A.l+B.l,A.u+B.u,A.rows,B.cols),
             A,B)
end


Base.convert(::Type{BlockBandedMatrix},B::Matrix) =
    BlockBandedMatrix(vec(B),0,0,[size(B,1)],[size(B,2)])

Base.convert(::Type{BlockBandedMatrix},B::BandedMatrix) =
    if isdiag(B)
        BlockBandedMatrix(copy(vec(B.data)),0,0,ones(Int,size(B,1)),ones(Int,size(B,2)))
    else
        BlockBandedMatrix(Matrix(B))
    end


## back substitution
#TODO: don't auto-pad
function trtrs!(::Type{Val{'U'}},A::BlockBandedMatrix,u::Vector)
    # When blocks are square, use LAPACK trtrs!
    mn=min(length(A.rows),length(A.cols))
    if A.rows[1:mn] == A.cols[1:mn]
        blockbanded_squareblocks_trtrs!(A,u)
    else
        blockbanded_rectblocks_trtrs!(A,u)
    end
end



function blockbanded_squareblocks_trtrs!(A::BlockBandedMatrix,u::Vector)
    if size(A,1) < size(u,1)
        throw(BoundsError())
    end
    n=size(u,1)
    N=Block(A.rowblocks[n])

    kr1=blockrows(A,N)
    b=n-kr1[1]+1
    kr1=kr1[1]:n

    trtrs!('U','N','N',view(A,N[1:b],N[1:b]),view(u,kr1))

    for K=N-1:-1:Block(1)
        kr=blockrows(A,K)
        for J=min(N,blockrowstop(A,K)):-1:K+1
            if J==N  # need to take into account zeros
                gemv!('N',-one(eltype(A)),view(A,K,N[1:b]),view(u,kr1),one(eltype(A)),view(u,kr))
            else
                gemv!('N',-one(eltype(A)),view(A,K,J),view(u,blockcols(A,J)),one(eltype(A)),view(u,kr))
            end
        end
        trtrs!('U','N','N',view(A,K,K),view(u,kr))
    end

    u
end

function blockbanded_rectblocks_trtrs!{T}(R::BlockBandedMatrix{T},b::Vector)
    n=n_end=length(b)
    K_diag=N=Block(R.rowblocks[n])
    J_diag=M=Block(R.colblocks[n])

    while n > 0
        B_diag = view(R,K_diag,J_diag)

        kr = blockrows(R,K_diag)
        jr = blockcols(R,J_diag)


        k = n-kr[1]+1
        j = n-jr[1]+1

        skr = max(1,k-j+1):k   # range in the sub block
        sjr = max(1,j-k+1):j   # range in the sub block

        kr2 = kr[skr]  # diagonal rows/cols we are working with

        for J = min(M,blockrowstop(R,K_diag)):-1:J_diag+1
            B=view(R,K_diag,J)
            Sjr = blockcols(R,J)

            if J==M
                Sjr = Sjr[1]:n_end  # The sub rows of the rhs we will multiply
                gemv!('N',-one(T),view(B,skr,1:length(Sjr)),
                                    view(b,Sjr),one(T),view(b,kr2))
            else  # can use all columns
                gemv!('N',-one(T),view(B,skr,:),
                                    view(b,Sjr),one(T),view(b,kr2))
            end
        end

        if J_diag ≠ M && sjr[end] ≠ size(B_diag,2)
            # subtract non-triangular columns
            sjr2 = sjr[end]+1:size(B_diag,2)
            gemv!('N',-one(T),view(B_diag,skr,sjr2),
                            view(b,sjr2 + jr[1]-1),one(T),view(b,kr2))
        elseif J_diag == M && sjr[end] ≠ size(B_diag,2)
            # subtract non-triangular columns
            Sjr = jr[1]+sjr[end]:n_end
            gemv!('N',-one(T),view(B_diag,skr,sjr[end]+1:sjr[end]+length(Sjr)),
                            view(b,Sjr),one(T),view(b,kr2))
        end

        trtrs!('U','N','N',view(B_diag,skr,sjr),view(b,kr2))

        if k == j
            K_diag -= 1
            J_diag -= 1
        elseif j < k
            J_diag -= 1
        else # if k < j
            K_diag -= 1
        end

        n = kr2[1]-1
    end
    b
end


function trtrs!{T}(A::BlockBandedMatrix{T},u::Matrix)
    if size(A,1) < size(u,1)
        throw(BoundsError())
    end
    n=size(u,1)
    N=Block(A.rowblocks[n])

    kr1=blockrows(A,N)
    b=n-kr1[1]+1
    kr1=kr1[1]:n

    trtrs!('U','N','N',view(A,N[1:b],N[1:b]),view(u,kr1,:))

    for K=N-1:-1:Block(1)
        kr=blockrows(A,K)
        for J=min(N,blockrowstop(A,K)):-1:K+1
            if J==N  # need to take into account zeros
                gemm!('N',-one(T),view(A,K,N[1:b]),view(u,kr1,:),one(T),view(u,kr,:))
            else
                gemm!('N',-one(T),view(A,K,J),view(u,blockcols(A,J),:),one(T),view(u,kr,:))
            end
        end
        trtrs!('U','N','N',view(A,K,K),view(u,kr,:))
    end

    u
end



## reindex

index_shape(A::AbstractArray,  I::SubBlock)    = index_shape(A,I.sub)


function reindex{T}(A::SubArray{T,2}, B::Tuple{UnitRange{Int},Any}, kj::Tuple{Block,Any})
    K = first(kj)
    cols = blockrows(parent(A),K)
    (SubBlock(K,(cols ∩ first(B)) - first(cols) + 1),reindex(A,tail(B),tail(kj))[1])
end

function reindex{T}(A::SubArray{T,2}, B::Tuple{UnitRange{Int}}, j::Tuple{Block})
    J = first(j)
    cols = blockcols(parent(A),J)
    (SubBlock(J,(cols ∩ first(B)) - first(cols) + 1),)
end

function reindex{T,SSB<:SubBlock}(A::SubArray{T,2}, B::Tuple{UnitRange{Int},Any}, kj::Tuple{SSB,Any})
    SB = reindex(A,B,(block(kj[1]),kj[2]))
    (SubBlock(block(SB[1]),SB[1].sub[kj[1].sub]),SB[2])
end
function reindex{T,SSB<:SubBlock}(A::SubArray{T,2}, B::Tuple{UnitRange{Int}}, j::Tuple{SSB})
    SB = reindex(A,B,(block(j[1]),))
    (SubBlock(block(SB[1]),SB[1].sub[kj[1].sub]),)
end


function reindex{T}(A::SubArray{T,2}, B::Tuple{Block,Any}, kj::Tuple{Int,Any})
    SB = first(B)
    k2 = first(kj)
    (blockrows(parent(A),block(SB))[1]+k2-1,reindex(A,tail(B),tail(kj))[1])
end
function reindex{T}(A::SubArray{T,2}, B::Tuple{Block}, j::Tuple{Int})
    SB = first(B)
    j2 = first(j)
    (blockcols(parent(A),block(SB))[1]+j2-1,)
end

function reindex{T,SSB<:SubBlock}(A::SubArray{T,2}, B::Tuple{SSB,Any}, kj::Tuple{Int,Any})
    SB = first(B)
    k2 = reindex(nothing,(SB.sub,),(first(kj),))[1]
    (blockrows(parent(A),block(SB))[1]+k2-1,reindex(A,tail(B),tail(kj))[1])
end
function reindex{T,SSB<:SubBlock}(A::SubArray{T,2}, B::Tuple{SSB}, j::Tuple{Int})
    SB = first(B)
    j2 = reindex(nothing,(SB.sub,),(first(j),))[1]
    (blockcols(parent(A),block(SB))[1]+j2-1,)
end


## view


if VERSION ≥ v"0.6-"
    Base.view(A::AbstractBlockBandedMatrix,K::Union{Block,SubBlock,UnitRange{Block}},J::Union{Block,SubBlock,UnitRange{Block}}) =
        SubArray(A, (K,J))
else
    Base.view(A::AbstractBlockBandedMatrix,K::Union{Block,SubBlock,UnitRange{Block}},J::Union{Block,SubBlock,UnitRange{Block}}) =
        SubArray(A, (K,J), subblocksize(A,K,J))
end

Base.view(A::SubBandedBlockSubBlock,::Colon,::Colon) =
    view(parent(A),parentindexes(A)[1],parentindexes(A)[2])
Base.view(A::SubBandedBlockSubBlock,::Colon,J::Union{AbstractArray,Real}) =
    view(parent(A),parentindexes(A)[1],parentindexes(A)[2][J])
Base.view(A::SubBandedBlockSubBlock,K::Union{AbstractArray,Real},::Colon) =
    view(parent(A),parentindexes(A)[1][K],parentindexes(A)[2])
Base.view(A::SubBandedBlockSubBlock,K::Union{AbstractArray,Real},J::Union{AbstractArray,Real}) =
    view(parent(A),parentindexes(A)[1][K],parentindexes(A)[2][J])
Base.view(A::SubBandedBlockSubBlock,K,J) =
    view(parent(A),parentindexes(A)[1][K],parentindexes(A)[2][J])

Base.view(A::SubBandedBlockRange,K::Union{Block,SubBlock,UnitRange{Block}},J::Union{Block,SubBlock,UnitRange{Block}}) =
    view(parent(A),parentindexes(A)[1][Int.(K)],parentindexes(A)[2][Int.(J)])

Base.view{T,BBM<:AbstractBlockBandedMatrix}(A::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}}},k::Block,j::Block) =
    view(parent(A),reindex(A,parentindexes(A),(k,j))...)
