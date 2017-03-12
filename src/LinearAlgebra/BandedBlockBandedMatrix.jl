# Represents a block banded matrix with banded blocks
#   similar to BandedMatrix{BandedMatrix{T}}
type BandedBlockBandedMatrix{T,RI,CI} <: AbstractBlockBandedMatrix{T,BandedMatrix{T}}
    data::Matrix{T}

    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    λ::Int  # sub lower bandwidth
    μ::Int  # sub upper bandwidth

    rows::RI   # the size of the row blocks
    cols::CI  # the size of the column blocks

    rowblocks::Vector{Int}
    colblocks::Vector{Int}

    function (::Type{BandedBlockBandedMatrix{T,RI,CI}}){T,RI,CI}(data::Matrix{T},l,u,λ,μ,rows,cols)
        @assert size(data,1) == λ+μ+1
        @assert size(data,2) == (l+u+1)*sum(cols)
        new{T,RI,CI}(data,l,u,λ,μ,rows,cols,blocklookup(rows),blocklookup(cols))
    end
end

BandedBlockBandedMatrix(data::Matrix,l,u,λ,μ,rows,cols) =
    BandedBlockBandedMatrix{eltype(data),typeof(rows),typeof(cols)}(data,l,u,λ,μ,rows,cols)

BandedBlockBandedMatrix{T}(::Type{T},l,u,λ,μ,rows,cols) =
    BandedBlockBandedMatrix(Matrix{T}(λ+μ+1,(l+u+1)*sum(cols)),l,u,λ,μ,rows,cols)

for FUNC in (:zeros,:rand,:ones)
    BFUNC = parse("bbb"*string(FUNC))
    @eval function $BFUNC{T}(::Type{T},l,u,λ,μ,rows,cols)
        data = $FUNC(T,λ+μ+1,(l+u+1)*sum(cols))::Matrix{T}
        BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols)
    end
end

function Base.convert(::Type{BandedBlockBandedMatrix},Y::AbstractMatrix)
    if !isbandedblockbanded(Y)
        error("$(typeof(Y)) is not banded block banded")
    end
    ret = bbbzeros(eltype(Y),blockbandwidth(Y,1),blockbandwidth(Y,2),
                        subblockbandwidth(Y,1),subblockbandwidth(Y,2),
                        rowblocklengths(Y),colblocklengths(Y))
    BLAS.axpy!(one(eltype(Y)),Y,ret)
end


function BandedMatrix(B::BandedBlockBandedMatrix)
    if length(B.rows) == length(B.cols) == 1
        BandedMatrix(view(B,Block(1),Block(1)))
    elseif all(x->x==1,B.rows) && all(x->x==1,B.cols)
        BandedMatrix(B.data,length(B.rows),B.l,B.u)
    else
        error("$B is not a banded matrix")
    end
end

isbandedblockbanded(::) = false
isbandedblockbanded(::BandedBlockBandedMatrix) = true


Base.isdiag(A::BandedBlockBandedMatrix) = A.λ == A.μ == A.l == A.u


zeroblock(X::BandedBlockBandedMatrix,K::Block,J::Block) = bzeros(eltype(X),X.rows[K.K],X.cols[J.K],X.λ,X.μ)


@compat const BandedBlockBandedBlock{T,U,V} = SubArray{T,2,BandedBlockBandedMatrix{T,U,V},Tuple{Block,Block},false}
@compat const BandedBlockBandedSubBlock{T,U,V} = SubArray{T,2,BandedBlockBandedMatrix{T,U,V},Tuple{SubBlock{UnitRange{Int}},SubBlock{UnitRange{Int}}},false}
@compat const SubBandedBlockBandedRange{T,BBM<:BandedBlockBandedMatrix} = SubArray{T,2,BBM,Tuple{UnitRange{Block},UnitRange{Block}},false}
@compat const BLASBandedMatrix2{T,A,I} = Union{BandedBlockBandedBlock{T,A,I},BandedMatrices.BLASBandedMatrix{T}}

isbandedblockbanded(::SubBandedBlockBandedRange) = true
isbandedblockbanded{T,BBM<:BandedBlockBandedMatrix}(::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}},false}) = true

BandedMatrices.isbanded{T,U,V}(::BandedBlockBandedBlock{T,U,V}) = true

function Base.view(A::BandedBlockBandedMatrix,K::Union{Block,SubBlock},J::Union{Block,SubBlock})
    m = A.cols[block(J).K]
    # we obuse the extra variables to store some data
    SubArray{eltype(A),2,typeof(A),Tuple{typeof(K),typeof(J)},false}(A, (K,J),
                sum(A.cols[1:block(J).K-1])*(A.l+A.u+1)+(block(K).K-block(J).K+A.u)*m+1,m)
end



function subblockbandinds{T,BBM<:BandedBlockBandedMatrix}(S::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}},false},k)
    P = parent(S)
    kr,jr = parentindexes(S)
    if k == 1
        JR1 = P.colblocks[jr[1]]
        j_st = jr[1]-blockcols(P,JR1)[1]+1
        λ = subblockbandwidth(P,1)
        -(λ + j_st - 1)
    else
        KR1 = P.rowblocks[kr[1]]
        # amount that a block is shift
        k_st = kr[1]-blockrows(P,KR1)[1]+1
        μ = subblockbandwidth(P,2)
        μ + k_st - 1
    end
end

# the bandwidths will grow from a truncated matrix








## Bandedmatrix support

@inline BandedMatrices.leadingdimension{T,U,V}(S::BandedBlockBandedBlock{T,U,V}) = stride(parent(S).data,2)
BandedMatrices.bandwidth{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, k::Int) = k==1 ? parent(S).λ : parent(S).μ


@inline BandedMatrices.leadingdimension{T,U,V}(S::BandedBlockBandedSubBlock{T,U,V}) = stride(parent(S).data,2)
function BandedMatrices.bandwidth{T,U,V}(S::BandedBlockBandedSubBlock{T,U,V}, k::Int)
    sh = first(parentindexes(S)[1].sub)-first(parentindexes(S)[2].sub)
    k==1 ? parent(S).λ-sh : parent(S).μ+sh
end


subblockbandwidths(K) = -subblockbandinds(K,1),subblockbandinds(K,2)
subblockbandinds(K) = subblockbandinds(K,1),subblockbandinds(K,2)
subblockbandwidth(K,k::Integer) = k==1?-subblockbandinds(K,k):subblockbandinds(K,k)
subblockbandinds(K::BandedBlockBandedMatrix,k::Integer) = k==1 ? -K.λ : K.μ





@inline function inbands_getindex{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    u=A.μ
    j_sh  = j+col-1  # shift column to get correct column
    A.data[u + k - j + 1, j_sh]
end

@inline function inbands_setindex!{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, v, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    u=A.μ
    j_sh  = j+col-1  # shift column to get correct column
    @inbounds A.data[u + k - j + 1, j_sh] = convert(T, v)::T
    v
end

@inline function inbands_getindex{T,U,V}(S::BandedBlockBandedSubBlock{T,U,V}, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    u=A.μ
    SK,SJ = parentindexes(S)
    k,j = SK.sub[k],SJ.sub[j]
    j_sh  = j+col-1  # shift column to get correct column
    A.data[u + k - j + 1, j_sh]
end

@inline function inbands_setindex!{T,U,V}(S::BandedBlockBandedSubBlock{T,U,V}, v, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    u=A.μ
    SK,SJ = parentindexes(S)
    k,j = SK.sub[k],SJ.sub[j]
    j_sh  = j+col-1  # shift column to get correct column
    @inbounds A.data[u + k - j + 1, j_sh] = convert(T, v)::T
    v
end


function getindex{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    BandedMatrices.banded_getindex(view(A.data,:,col:col+S.stride1-1),A.λ,A.μ,k,j)
end

function setindex!{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, v, k::Int, j::Int)
    A = parent(S)
    col = S.offset1 # first col of current block
    BandedMatrices.banded_setindex!(view(A.data,:,col:col+S.stride1-1),A.λ,A.μ,v,k,j)
end


getindex(A::BandedBlockBandedMatrix,kr::UnitRange{Int},jr::UnitRange{Int}) =
    BandedBlockBandedMatrix(view(A,kr,jr))


## algebra



BLAS.axpy!{T,U,V,T2,U2,V2}(α,A::BandedBlockBandedBlock{T,U,V},B::BandedBlockBandedBlock{T2,U2,V2}) =
    BandedMatrices.banded_axpy!(α,A,B)

BLAS.axpy!{T,U,V,T2,U2,V2}(α,A::BandedBlockBandedSubBlock{T,U,V},B::BandedBlockBandedBlock{T2,U2,V2}) =
    BandedMatrices.banded_axpy!(α,A,B)

BLAS.axpy!{T,U,V,T2,U2,V2}(α,A::BandedBlockBandedSubBlock{T,U,V},B::BandedBlockBandedSubBlock{T2,U2,V2}) =
    BandedMatrices.banded_axpy!(α,A,B)


BLAS.axpy!{T,U,V}(α,A::BandedBlockBandedBlock{T,U,V},B::SubBandedBlockSubBlock) =
    BandedMatrices.banded_dense_axpy!(α,A,blockview(B))

BLAS.axpy!{T,U,V}(α,A::BandedBlockBandedBlock{T,U,V},B::AbstractMatrix) =
    BandedMatrices.banded_dense_axpy!(α,A,B)

Base.BLAS.axpy!(α,A::SubBandedBlockBandedRange,Y::SubBandedBlockBandedRange) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::SubBandedBlockBandedRange,Y::AbstractBlockMatrix) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::SubBandedBlockBandedRange) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::SubBandedBlockBandedRange,Y::SubBandedBlockRange) = block_axpy!(α,A,Y)
Base.BLAS.axpy!(α,A::SubBandedBlockRange,Y::SubBandedBlockBandedRange) = block_axpy!(α,A,Y)


## Convert routines

function Base.convert{T,U,V}(::Type{BandedMatrix{T}},S::BandedBlockBandedBlock{T,U,V})
    A = parent(S)
    col = S.offset1 # first col of current block
    BandedMatrix(A.data[:,col:col+S.stride1-1],A.rows[K.K],A.λ,A.μ)
end

function Base.pointer{T<:BlasFloat,U,V}(S::BandedBlockBandedBlock{T,U,V})
    A = parent(S)
    K,J = parentindexes(S)
    if K.K < J.K-A.u || K.K > J.K+A.l
        error("Cannot view zero blocks")
    end
    # column block K-J+A.u+1,J
    p=pointer(A.data)
    st=stride(A.data,2)
    sz=sizeof(T)
    col=S.offset1
    p+(col-1)*st*sz
end

leadingdimension{T<:BlasFloat,U,V}(B::BandedBlockBandedSubBlock{T,U,V}) =
    stride(parent(B).data,2)

function Base.pointer{T<:BlasFloat,U,V}(B::BandedBlockBandedSubBlock{T,U,V})
    p = pointer(parentblock(B))
    p+leadingdimension(B)*(first(parentindexes(B)[2].sub)-1)*sizeof(T)
end


*{T,U,V}(A::BandedBlockBandedBlock{T,U,V},B::BandedBlockBandedBlock{T,U,V}) = BandedMatrices.banded_A_mul_B(A,b)
*{T,U,V}(A::BandedBlockBandedBlock{T,U,V},b::AbstractVector{T}) =
    BandedMatrices.banded_A_mul_B!(Vector{T}(size(A,1)),A,b)


Base.A_mul_B!{T,U,V}(c::AbstractVector,A::BandedBlockBandedBlock{T,U,V},b::AbstractVector) =
    banded_A_mul_B!(c,A,b)

αA_mul_B_plus_βC!{T,U,V}(α,A::BandedBlockBandedBlock{T,U,V},x::AbstractVector,β,y::AbstractVector) =
    BandedMatrices.gbmv!('N',α,A,x,β,y)

αA_mul_B_plus_βC!{T,U,V}(α,A::BandedBlockBandedSubBlock{T,U,V},x::AbstractVector,β,y::AbstractVector) =
    BandedMatrices.gbmv!('N',α,A,x,β,y)


αA_mul_B_plus_βC!{T,U,V}(α,A::BLASBandedMatrix2{T,U,V},B::BLASBandedMatrix2{T,U,V},β,C::BLASBandedMatrix2{T,U,V}) =
    BandedMatrices.gbmm!(α,A,B,β,C)


αA_mul_B_plus_βC!{T,BBM<:BandedBlockBandedMatrix}(α,A::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}},false},
                                                  x::AbstractVector,β,y::AbstractVector) =
    blockmatrix_αA_mul_B_plus_βC!(α,A,x,β,y)

Base.A_mul_B!{T,BBM<:BandedBlockBandedMatrix}(y::Vector,A::SubArray{T,2,BBM,Tuple{UnitRange{Int},UnitRange{Int}},false},b::Vector) =
    αA_mul_B_plus_βC!(one(eltype(A)),A,b,zero(eltype(y)),y)

## Algebra



function *{T<:Number,V<:Number}(A::BandedBlockBandedMatrix{T},
                                B::BandedBlockBandedMatrix{V})
    if A.cols ≠ B.rows
        # diagonal matrices can be converted
        if isdiag(B) && size(A,2) == size(B,1) == size(B,2)
            B = BandedBlockBandedMatrix(B.data,0,0,0,0,A.cols,A.cols)
        elseif isdiag(A) && size(A,2) == size(B,1) == size(A,1)
            A = BandedBlockBandedMatrix(A.data,0,0,0,0,B.rows,B.rows)
        else
            throw(DimensionMismatch("*"))
        end
    end
    n,m=size(A,1),size(B,2)

    A_mul_B!(BandedBlockBandedMatrix(promote_type(T,V),A.l+B.l,A.u+B.u,
                                     A.λ+B.λ,A.μ+B.μ,A.rows,B.cols),
             A,B)
end



# convert




Base.convert(::Type{BandedBlockBandedMatrix},B::BandedMatrix) =
    if isdiag(B)
        BandedBlockBandedMatrix(copy(B.data),0,0,0,0,ones(Int,size(B,1)),ones(Int,size(B,2)))
    else
        BandedBlockBandedMatrix(copy(B.data),0,0,B.l,B.u,[size(B,1)],[size(B,2)])
    end
