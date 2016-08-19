# Represents a block banded matrix with banded blocks


function blocklookup(rows)
    rowblocks=Array(Int,0)
    for ν in eachindex(rows), k in 1:rows[ν]
        push!(rowblocks,ν)
    end
    rowblocks
end



immutable BandedBlockBandedMatrix{T,RI,CI} <: AbstractMatrix{T}
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


Base.size(A::BandedBlockBandedMatrix) = sum(A.rows),sum(A.cols)

Base.isdiag(A::BandedBlockBandedMatrix) = A.λ == A.μ == A.l == A.u

blocksize(A::BandedBlockBandedMatrix) = length(A.rows),length(A.cols)
blocksize(A::BandedBlockBandedMatrix,k::Int) = k==1?length(A.rows):length(A.cols)


blockcolrange(A::BandedBlockBandedMatrix,J::Int) = max(1,J-A.u):min(length(A.rows),J+A.l)


function getblock(A::BandedBlockBandedMatrix,K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        zeros(eltype(A),A.rows[K],A.cols[J])
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(l+u+1)  # number of columns before current block
        cols=S+(K-J+A.u)A.cols[J]+1:S+(K-J+A.u+1)A.cols[J]
        BandedMatrix(A.data[:,cols],
                        A.rows[K],A.λ,A.μ)
    end
end

function viewblock{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        error("Cannot view zero blocks")
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(A.l+A.u+1)  # number of columns before current block
        p=pointer(A.data)
        st=stride(A.data,2)
        sz=sizeof(p)
        cols=S+(K-J+A.u)A.cols[J]+1:S+(K-J+A.u+1)A.cols[J]

        p+=(first(cols)-1)*st*sz
        BandedMatrix(unsafe_wrap(Array,p,(st,length(cols))),
                        A.rows[K],A.λ,A.μ)
    end
end


function Base.getindex{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        zero(eltype(A))
    else
        k2=k-sum(A.rows[1:K-1])
        j2=j-sum(A.cols[1:J-1])
        viewblock(A,K,J)[k2,j2]
    end
end

function Base.setindex!(A::BandedBlockBandedMatrix,v,k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    k2=k-sum(A.rows[1:K-1])
    j2=j-sum(A.cols[1:J-1])
    setindex!(viewblock(A,K,J),v,k2,j2)
end

Base.linearindexing{BBBM<:BandedBlockBandedMatrix}(::Type{BBBM}) =
    Base.LinearSlow()


function gbmv!(α,A::BandedBlockBandedMatrix,x::Vector,β,y::Vector)
    if length(x) != size(A,2) || length(y) != size(A,1)
        throw(BoundsError())
    end

    BLAS.scal!(length(y),β,y,1)
    o=one(eltype(y))

    for J=1:blocksize(A,2)
        jr=blockcols(A,J)
        for K=blockcolrange(A,J)
            kr=blockrows(A,K)
            B=viewblock(A,K,J)
            gbmv!('N',α,B,view(x,jr),o,view(y,kr))
        end
    end
    y
end

Base.A_mul_B!(y::Vector,A::BandedBlockBandedMatrix,b::Vector) =
    gbmv!(one(eltype(A)),A,b,zero(eltype(y)),y)


blockrows(A::BandedBlockBandedMatrix,K::Int) =
    sum(A.rows[1:K-1]) + (1:A.rows[K])
blockcols(A::BandedBlockBandedMatrix,J::Int) =
    sum(A.cols[1:J-1]) + (1:A.cols[J])

function Base.BLAS.axpy!(α,A::BandedBlockBandedMatrix,Y::Matrix)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=1:blocksize(A,2)
        jr=blockcols(A,J)
        for K=blockcolrange(A,J)
            kr=blockrows(A,K)
            BLAS.axpy!(α,viewblock(A,K,J),view(Y,kr,jr))
        end
    end
    Y
end

function Base.BLAS.axpy!(α,A::BandedBlockBandedMatrix,Y::BandedBlockBandedMatrix)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=1:blocksize(A,2), K=blockcolrange(A,J)
        BLAS.axpy!(α,viewblock(A,K,J),viewblock(Y,K,J))
    end
    Y
end

Base.convert(::Type{Matrix},A::BandedBlockBandedMatrix) =
    BLAS.axpy!(one(eltype(A)),A,zeros(eltype(A),size(A,1),size(A,2)))

Base.full(S::BandedBlockBandedMatrix) = convert(Matrix, S)



## Algebra


function Base.A_mul_B!(Y::BandedBlockBandedMatrix,
                       A::BandedBlockBandedMatrix,
                       B::BandedBlockBandedMatrix)
    T=eltype(Y)
    BLAS.scal!(length(Y.data),zero(T),Y.data,1)
    o=one(T)
    for J=1:blocksize(B,2),N=blockcolrange(B,J),K=blockcolrange(A,N)
        gbmm!(o,viewblock(A,K,N),viewblock(B,N,J),o,viewblock(Y,K,J))
    end
    Y
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
