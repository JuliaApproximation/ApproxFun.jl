# Represents a block banded matrix with banded blocks
#   similar to BandedMatrix{BandedMatrix{T}}
type BandedBlockBandedMatrix{T,RI,CI} <: AbstractBandedBlockMatrix{T}
    data::Matrix{T}

    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    λ::Int  # sub lower bandwidth
    μ::Int  # sub upper bandwidth

    rows::RI   # the size of the row blocks
    cols::CI  # the size of the column blocks

    rowblocks::Vector{Int}
    colblocks::Vector{Int}

    function BandedBlockBandedMatrix(data::Matrix{T},l,u,λ,μ,rows,cols)
        @assert size(data,1) == λ+μ+1
        @assert size(data,2) == (l+u+1)*sum(cols)
        new(data,l,u,λ,μ,rows,cols,blocklookup(rows),blocklookup(cols))
    end
end

BandedBlockBandedMatrix(data::Matrix,l,u,λ,μ,rows,cols) =
    BandedBlockBandedMatrix{eltype(data),typeof(rows),typeof(cols)}(data,l,u,λ,μ,rows,cols)

BandedBlockBandedMatrix{T}(::Type{T},l,u,λ,μ,rows,cols) =
    BandedBlockBandedMatrix(Array(T,λ+μ+1,(l+u+1)*sum(cols)),l,u,λ,μ,rows,cols)

for FUNC in (:zeros,:rand,:ones)
    BFUNC = parse("bbb"*string(FUNC))
    @eval $BFUNC{T}(::Type{T},l,u,λ,μ,rows,cols) =
        BandedBlockBandedMatrix($FUNC(T,λ+μ+1,(l+u+1)*sum(cols)),l,u,λ,μ,rows,cols)
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


Base.isdiag(A::BandedBlockBandedMatrix) = A.λ == A.μ == A.l == A.u



typealias BandedBlockBandedBlock{T,U,V} SubArray{T,2,BandedBlockBandedMatrix{T,U,V},Tuple{Block,Block},false}
typealias BLASBandedMatrix2{T,A,I} Union{BandedBlockBandedBlock{T,A,I},BandedMatrices.BLASBandedMatrix{T}}


Base.view(A::BandedBlockBandedMatrix,K::Block,J::Block) =
    SubArray(A, (K,J), (A.rows[K.K],A.cols[J.K]))
Base.indices{T,U,V}(S::BandedBlockBandedBlock{T,U,V}) =
    (Base.OneTo(parent(S).rows[parentindexes(S)[1].K]),
     Base.OneTo(parent(S).cols[parentindexes(S)[2].K]))



@inline BandedMatrices.leadingdimension{T,U,V}(S::BandedBlockBandedBlock{T,U,V}) = stride(parent(S).data,2)
BandedMatrices.bandwidth{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, k::Int) = k==1 ? parent(S).λ : parent(S).μ

function block_pointer{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        error("Cannot view zero blocks")
    end
    # column block K-J+A.u+1,J
    S=sum(A.cols[1:J-1])*(A.l+A.u+1)  # number of columns before current block
    p=pointer(A.data)
    st=stride(A.data,2)
    sz=sizeof(T)
    col=S+(K-J+A.u)A.cols[J]+1
    p+(col-1)*st*sz
end


function getindex{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, k::Int, j::Int)
    A = parent(S)
    K,J = parentindexes(S)
    m = A.cols[J.K]

    S = sum(A.cols[1:J.K-1])*(A.l+A.u+1)  # number of columns before current block
    col = S+(K.K-J.K+A.u)*m+1

    BandedMatrices.banded_getindex(view(A.data,:,col:col+m-1),A.λ,A.μ,k,j)
end

function setindex!{T,U,V}(S::BandedBlockBandedBlock{T,U,V}, v, k::Int, j::Int)
    A = parent(S)
    K,J = parentindexes(S)
    m = A.cols[J.K]

    S = sum(A.cols[1:J.K-1])*(A.l+A.u+1)  # number of columns before current block
    col = S+(K.K-J.K+A.u)*m+1

    BandedMatrices.banded_setindex!(view(A.data,:,col:col+m-1),A.λ,A.μ,v,k,j)
end

function Base.convert{T,U,V}(::Type{BandedMatrix{T}},S::BandedBlockBandedBlock{T,U,V})
    A = parent(S)
    K,J = parentindexes(S)
    m = A.cols[J.K]

    S = sum(A.cols[1:J.K-1])*(A.l+A.u+1)  # number of columns before current block
    col = S+(K.K-J.K+A.u)*m+1
    BandedMatrix(A.data[:,col:col+m-1],A.rows[K.K],A.λ,A.μ)
end

function Base.pointer{T<:BlasFloat,U,V}(S::BandedBlockBandedBlock{T,U,V})
    A = parent(S)
    K,J = parentindexes(S)
    block_pointer(A,K.K,J.K)
end


*{T,U,V}(A::BandedBlockBandedBlock{T,U,V},B::BandedBlockBandedBlock{T,U,V}) = BandedMatrices.banded_A_mul_B(A,b)
*{T,U,V}(A::BandedBlockBandedBlock{T,U,V},b::AbstractVector{T}) = BandedMatrices.banded_A_mul_B!(Array(T,length(b)),A,b)


Base.A_mul_B!{T,U,V}(c::AbstractVector,A::BandedBlockBandedBlock{T,U,V},b::AbstractVector) =
    banded_A_mul_B!(c,A,b)

αA_mul_B_plus_βC!{T,U,V}(α,A::BandedBlockBandedBlock{T,U,V},x::AbstractVector,β,y::AbstractVector) =
    BandedMatrices.gbmv!('N',α,A,x,β,y)
αA_mul_B_plus_βC!{T,U,V}(α,A::BLASBandedMatrix2{T,U,V},B::BLASBandedMatrix2{T,U,V},β,C::BLASBandedMatrix2{T,U,V}) =
    BandedMatrices.gbmm!(α,A,B,β,C)


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



Base.convert(::Type{BandedBlockBandedMatrix},B::BandedMatrix) =
    if isdiag(B)
        BandedBlockBandedMatrix(copy(B.data),0,0,0,0,ones(Int,size(B,1)),ones(Int,size(B,2)))
    else
        BandedBlockBandedMatrix(copy(B.data),0,0,B.l,B.u,[size(B,1)],[size(B,2)])
    end
