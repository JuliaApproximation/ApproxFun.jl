# Represents a block banded matrix with banded blocks
#   similar to BandedMatrix{BandedMatrix{T}}
immutable BandedBlockBandedMatrix{T,RI,CI} <: AbstractBandedBlockMatrix{T}
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



Base.isdiag(A::BandedBlockBandedMatrix) = A.λ == A.μ == A.l == A.u


function viewblock{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        error("Cannot view zero blocks")
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(A.l+A.u+1)  # number of columns before current block
        p=pointer(A.data)
        st=stride(A.data,2)
        sz=sizeof(T)
        cols=S+(K-J+A.u)A.cols[J]+1:S+(K-J+A.u+1)A.cols[J]

        p+=(first(cols)-1)*st*sz
        BandedMatrix(unsafe_wrap(Array,p,(st,length(cols))),
                        A.rows[K],A.λ,A.μ)
    end
end


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
