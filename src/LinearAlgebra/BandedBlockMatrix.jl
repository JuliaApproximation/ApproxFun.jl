# Takes in a list of lengths, and allows lookup of which block an entry is in
function blocklookup(rows)
    rowblocks=Array(Int,0)
    for ν in eachindex(rows), k in 1:rows[ν]
        push!(rowblocks,ν)
    end
    rowblocks
end


abstract AbstractBlockMatrix{T} <: AbstractMatrix{T}
abstract AbstractBandedBlockMatrix{T} <: AbstractBlockMatrix{T}

Base.size(A::AbstractBlockMatrix) = sum(A.rows),sum(A.cols)


blocksize(A::AbstractBlockMatrix) = length(A.rows),length(A.cols)
blocksize(A::AbstractBlockMatrix,k::Int) = k==1?length(A.rows):length(A.cols)

# these give the block rows corresponding to the J-th column block
blockcolstart(A::AbstractBandedBlockMatrix,J::Int) = max(1,J-A.u)
blockcolstop(A::AbstractBandedBlockMatrix,J::Int) = min(length(A.rows),J+A.l)
blockcolrange(A::AbstractBandedBlockMatrix,J::Int) = blockcolstart(A,J):blockcolstop(A,J)

blockrowstart(A::AbstractBandedBlockMatrix,K::Int) = max(1,K-A.l)
blockrowstop(A::AbstractBandedBlockMatrix,K::Int) = min(length(A.cols),K+A.u)
blockrowrange(A::AbstractBandedBlockMatrix,K::Int) = blockrowstart(A,K):blockrowstop(A,K)


# give the rows/columns of a block, as a range
blockrows(A::AbstractBlockMatrix,K::Int) = sum(A.rows[1:K-1]) + (1:A.rows[K])
blockcols(A::AbstractBlockMatrix,J::Int) = sum(A.cols[1:J-1]) + (1:A.cols[J])

colstop(A::AbstractBandedBlockMatrix,k::Int) = sum(A.rows[1:min(A.colblocks[k]+A.l,length(A.rows))])
colstart(A::AbstractBandedBlockMatrix,k::Int) = sum(A.rows[1:A.colblocks[k]-A.u-1])+1

rowstop(A::AbstractBandedBlockMatrix,k::Int) = sum(A.cols[1:min(A.rowblocks[k]+A.u,length(A.cols))])
rowstart(A::AbstractBandedBlockMatrix,k::Int) = sum(A.cols[1:A.rowblocks[k]-A.l-1])+1

Base.convert(::Type{Matrix},A::AbstractBlockMatrix) =
    BLAS.axpy!(one(eltype(A)),A,zeros(eltype(A),size(A,1),size(A,2)))

Base.full(S::AbstractBlockMatrix) = convert(Matrix, S)

function getblock(A::AbstractBandedBlockMatrix,K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        zeros(eltype(A),A.rows[K],A.cols[J])
    else
        # column block K-J+A.u+1,J
        deepcopy(viewblock(A,K,J))
    end
end


function Base.getindex(A::AbstractBandedBlockMatrix,k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        zero(eltype(A))
    else
        k2=k-sum(A.rows[1:K-1])
        j2=j-sum(A.cols[1:J-1])
        viewblock(A,K,J)[k2,j2]
    end
end

function Base.setindex!(A::AbstractBandedBlockMatrix,v,k::Int,j::Int)
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
        setindex!(viewblock(A,K,J),v,k2,j2)
    end
end

Base.linearindexing{BBBM<:AbstractBlockMatrix}(::Type{BBBM}) =
    Base.LinearSlow()


function αA_mul_B_plus_βC!(α,A::AbstractBlockMatrix,x::Vector,β,y::Vector)
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
            αA_mul_B_plus_βC!(α,B,view(x,jr),o,view(y,kr))
        end
    end
    y
end

Base.A_mul_B!(y::Vector,A::AbstractBlockMatrix,b::Vector) =
    αA_mul_B_plus_βC!(one(eltype(A)),A,b,zero(eltype(y)),y)


function Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::AbstractMatrix)
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

function Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::AbstractBlockMatrix)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=1:blocksize(A,2), K=blockcolrange(A,J)
        BLAS.axpy!(α,viewblock(A,K,J),viewblock(Y,K,J))
    end
    Y
end




## Algebra


function Base.A_mul_B!(Y::AbstractBlockMatrix,A::AbstractBlockMatrix,B::AbstractBlockMatrix)
    T=eltype(Y)
    BLAS.scal!(length(Y.data),zero(T),Y.data,1)
    o=one(T)
    for J=1:blocksize(B,2),N=blockcolrange(B,J),K=blockcolrange(A,N)
        αA_mul_B_plus_βC!(o,viewblock(A,K,N),viewblock(B,N,J),o,viewblock(Y,K,J))
    end
    Y
end


## BandedBlockMatrix routines


function bbm_blockstarts(rows,cols,l,u)
    numentries = 0
    N=length(rows)
    bs=BandedMatrix(Int,length(rows),length(cols),l,u)
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
type BandedBlockMatrix{T,RI,CI} <: AbstractBandedBlockMatrix{T}
    data::Vector{T}

    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    rows::RI   # the size of the row blocks
    cols::CI  # the size of the column blocks

    rowblocks::Vector{Int}
    colblocks::Vector{Int}

    blockstart::BandedMatrix{Int}  # gives shift of first entry of data for each block

    function BandedBlockMatrix(data::Vector{T},l,u,rows,cols)
        new(data,l,u,rows,cols,blocklookup(rows),blocklookup(cols),bbm_blockstarts(rows,cols,l,u))
    end
end






BandedBlockMatrix(data::Vector,l,u,rows,cols) =
    BandedBlockMatrix{eltype(data),typeof(rows),typeof(cols)}(data,l,u,rows,cols)

BandedBlockMatrix{T}(::Type{T},l,u,rows,cols) =
    BandedBlockMatrix(Array(T,bbm_numentries(rows,cols,l,u)),l,u,rows,cols)

for FUNC in (:zeros,:rand,:ones)
    BFUNC = parse("bb"*string(FUNC))
    @eval $BFUNC{T}(::Type{T},l,u,rows,cols) =
        BandedBlockMatrix($FUNC(T,bbm_numentries(rows,cols,l,u)),l,u,rows,cols)
end


#TODO: in 0.5, use view
function viewblock{T<:BlasFloat}(A::BandedBlockMatrix{T},K::Int,J::Int)
    if K < J-A.u || K > J+A.l
        error("Cannot view zero blocks")
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(A.l+A.u+1)  # number of columns before current block
        p=pointer(A.data)
        sz=sizeof(p)


        p+=A.blockstart[K,J]*sz
        unsafe_wrap(Array,p,(A.rows[K],A.cols[J]))
    end
end


function *{T<:Number,V<:Number}(A::BandedBlockMatrix{T},
                                B::BandedBlockMatrix{V})
    if A.cols ≠ B.rows
        throw(DimensionMismatch("*"))
    end
    n,m=size(A,1),size(B,2)

    A_mul_B!(BandedBlockMatrix(promote_type(T,V),A.l+B.l,A.u+B.u,A.rows,B.cols),
             A,B)
end


Base.convert(::Type{BandedBlockMatrix},B::Matrix) =
    BandedBlockMatrix(vec(B),0,0,[size(B,1)],[size(B,2)])

Base.convert(::Type{BandedBlockMatrix},B::BandedMatrix) =
    if isdiag(B)
        BandedBlockMatrix(copy(B.data),0,0,ones(Int,size(B,1)),ones(Int,size(B,2)))
    else
        BandedBlockMatrix(Matrix(B))
    end


## back substitution
#TODO: don't auto-pad
function trtrs!(::Type{Val{'U'}},A::BandedBlockMatrix,u::Vector)
    # When blocks are square, use LAPACK trtrs!
    mn=min(length(A.rows),length(A.cols))
    if A.rows[1:mn] == A.cols[1:mn]
        bandedblock_squareblocks_trtrs!(A,u)
    else
        bandedblock_rectblocks_trtrs!(A,u)
    end
end



function bandedblock_squareblocks_trtrs!(A::BandedBlockMatrix,u::Vector)
    if size(A,1) < size(u,1)
        throw(BoundsError())
    end
    n=size(u,1)
    N=A.rowblocks[n]

    kr1=blockrows(A,N)
    b=n-kr1[1]+1
    kr1=kr1[1]:n

    LAPACK.trtrs!('U','N','N',view(viewblock(A,N,N),1:b,1:b),view(u,kr1))

    for K=N-1:-1:1
        kr=blockrows(A,K)
        for J=min(N,blockrowstop(A,K)):-1:K+1
            if J==N  # need to take into account zeros
                BLAS.gemv!('N',-1.0,view(viewblock(A,K,N),:,1:b),view(u,kr1),1.0,view(u,kr))
            else
                BLAS.gemv!('N',-1.0,viewblock(A,K,J),view(u,blockcols(A,J)),1.0,view(u,kr))
            end
        end
        LAPACK.trtrs!('U','N','N',viewblock(A,K,K),view(u,kr))
    end

    u
end

function bandedblock_rectblocks_trtrs!(R::BandedBlockMatrix,b::Vector)
    n=n_end=length(b)
    K_diag=N=R.rowblocks[n]
    J_diag=M=R.colblocks[n]

    while n > 0
        B_diag = viewblock(R,K_diag,J_diag)

        kr = blockrows(R,K_diag)
        jr = blockcols(R,J_diag)


        k = n-kr[1]+1
        j = n-jr[1]+1

        skr = max(1,k-j+1):k   # range in the sub block
        sjr = max(1,j-k+1):j   # range in the sub block

        kr2 = kr[skr]  # diagonal rows/cols we are working with

        for J = min(M,blockrowstop(R,K_diag)):-1:J_diag+1
            B=viewblock(R,K_diag,J)
            Sjr = blockcols(R,J)

            if J==M
                Sjr = Sjr[1]:n_end  # The sub rows of the rhs we will multiply
                BLAS.gemv!('N',-1.0,view(B,skr,1:length(Sjr)),
                                    view(b,Sjr),1.0,view(b,kr2))
            else  # can use all columns
                BLAS.gemv!('N',-1.0,view(B,skr,:),
                                    view(b,Sjr),1.0,view(b,kr2))
            end
        end

        if J_diag ≠ M && sjr[end] ≠ size(B_diag,2)
            # subtract non-triangular columns
            sjr2 = sjr[end]+1:size(B_diag,2)
            BLAS.gemv!('N',-1.0,view(B_diag,skr,sjr2),
                            view(b,sjr2 + jr[1]-1),1.0,view(b,kr2))
        elseif J_diag == M && sjr[end] ≠ size(B_diag,2)
            # subtract non-triangular columns
            Sjr = jr[1]+sjr[end]:n_end
            BLAS.gemv!('N',-1.0,view(B_diag,skr,sjr[end]+1:sjr[end]+length(Sjr)),
                            view(b,Sjr),1.0,view(b,kr2))
        end

        LAPACK.trtrs!('U','N','N',view(B_diag,skr,sjr),view(b,kr2))

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


function trtrs!(A::BandedBlockMatrix,u::Matrix)
    if size(A,1) < size(u,1)
        throw(BoundsError())
    end
    n=size(u,1)
    N=A.rowblocks[n]

    kr1=blockrows(A,N)
    b=n-kr1[1]+1
    kr1=kr1[1]:n

    LAPACK.trtrs!('U','N','N',view(viewblock(A,N,N),1:b,1:b),view(u,kr1,:))

    for K=N-1:-1:1
        kr=blockrows(A,K)
        for J=min(N,blockrowstop(A,K)):-1:K+1
            if J==N  # need to take into account zeros
                BLAS.gemm!('N',-1.0,view(viewblock(A,K,N),:,1:b),view(u,kr1,:),1.0,view(u,kr,:))
            else
                BLAS.gemm!('N',-1.0,viewblock(A,K,J),view(u,blockcols(A,J),:),1.0,view(u,kr,:))
            end
        end
        LAPACK.trtrs!('U','N','N',viewblock(A,K,K),view(u,kr,:))
    end

    u
end
