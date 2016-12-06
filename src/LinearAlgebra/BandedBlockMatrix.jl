doc"""
`Block` is used for get index of a block matrix.  For example,
```julia
A[Block(1),Block(2)]
```
retrieves the 1 x 2 block
"""
immutable Block
    K::Int
end

block(B::Block) = B

Base.convert{T<:Integer}(::Type{T},B::Block) = convert(T,B.K)::T
Base.convert(::Type{Block},K::Integer) = Block(K)

for OP in (:(Base.one),:(Base.zero),:(-))
    @eval $OP(B::Block) = Block($OP(B.K))
end

for OP in (:(==),:(!=),:(Base.isless))
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

done(1,true)

# broadcasting
Base.broadcast(Int,r::StepRange{Block,Block}) = Int(r.start):Int(r.step):Int(r.stop)
Base.broadcast(Int,B::Block) = Int(B)


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
    rowblocks=Array(Int,0)
    for ν in eachindex(rows), k in 1:rows[ν]
        push!(rowblocks,ν)
    end
    rowblocks
end


abstract AbstractBlockMatrix{T} <: AbstractMatrix{T}
abstract AbstractBandedBlockMatrix{T} <: AbstractBlockMatrix{T}


getindex(A::AbstractBlockMatrix,K::Block,J::Block) = copy(view(A,K,J))
getindex(A::AbstractBlockMatrix,K::Block,j) = A[blockrows(A,K),j]
getindex(A::AbstractBlockMatrix,k,J::Block) = A[k,blockcols(A,J)]

Base.size(A::AbstractBlockMatrix) = sum(A.rows),sum(A.cols)


blocksize(A::AbstractBlockMatrix) = length(A.rows),length(A.cols)
blocksize(A::AbstractBlockMatrix,k::Int) = k==1?length(A.rows):length(A.cols)

# these give the block rows corresponding to the J-th column block
blockcolstart(A::AbstractBandedBlockMatrix,J::Int) = Block(max(1,J-A.u))
blockcolstop(A::AbstractBandedBlockMatrix,J::Int) = Block(min(length(A.rows),J+A.l))
blockcolrange(A::AbstractBandedBlockMatrix,J) = blockcolstart(A,J):blockcolstop(A,J)

blockrowstart(A::AbstractBandedBlockMatrix,K::Int) = Block(max(1,K-A.l))
blockrowstop(A::AbstractBandedBlockMatrix,K::Int) = Block(min(length(A.cols),K+A.u))
blockrowrange(A::AbstractBandedBlockMatrix,K) = blockrowstart(A,K):blockrowstop(A,K)

# give the rows/columns of a block, as a range
blockrows(A::AbstractBlockMatrix,K::Int) = sum(A.rows[1:K-1]) + (1:A.rows[K])
blockcols(A::AbstractBlockMatrix,J::Int) = sum(A.cols[1:J-1]) + (1:A.cols[J])

for Op in (:blockcolstart,:blockcolstop,:blockrowstart,:blockrowstop,:blockrows,:blockcols)
    @eval $Op(A::AbstractBlockMatrix,K::Block) = $Op(A,K.K)
end

colstop(A::AbstractBandedBlockMatrix,k::Int) = sum(A.rows[1:min(A.colblocks[k]+A.l,length(A.rows))])
function colstart(A::AbstractBandedBlockMatrix,k::Int)
    lr = A.colblocks[k]-A.u-1
    lr ≤ length(A.rows) ? sum(A.rows[1:lr])+1 : size(A,1)+1  # no columns
end

rowstop(A::AbstractBandedBlockMatrix,k::Int) = sum(A.cols[1:min(A.rowblocks[k]+A.u,length(A.cols))])
rowstart(A::AbstractBandedBlockMatrix,k::Int) = sum(A.cols[1:A.rowblocks[k]-A.l-1])+1

Base.convert(::Type{Matrix},A::AbstractBlockMatrix) =
    BLAS.axpy!(one(eltype(A)),A,zeros(eltype(A),size(A,1),size(A,2)))

Base.full(S::AbstractBlockMatrix) = convert(Matrix, S)


function Base.getindex(A::AbstractBandedBlockMatrix,k::Int,j::Int)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        zero(eltype(A))
    else
        k2=k-sum(A.rows[1:K-1])
        j2=j-sum(A.cols[1:J-1])
        view(A,Block(K),Block(J))[k2,j2]
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
        view(A,Block(K),Block(J))[k2,j2] = v
    end
end

Base.linearindexing{BBBM<:AbstractBlockMatrix}(::Type{BBBM}) =
    Base.LinearSlow()


function αA_mul_B_plus_βC!(α,A::AbstractBlockMatrix,x::AbstractVector,β,y::AbstractVector)
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
            BLAS.axpy!(α,view(A,Block(K),Block(J)),view(Y,kr,jr))
        end
    end
    Y
end

function Base.BLAS.axpy!(α,A::AbstractBlockMatrix,Y::AbstractBlockMatrix)
    if size(A) ≠ size(Y)
        throw(BoundsError())
    end

    for J=1:blocksize(A,2), K=blockcolrange(A,J)
        BLAS.axpy!(α,view(A,Block(K),Block(J)),view(Y,Block(K),Block(J)))
    end
    Y
end




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

    rowblocks::Vector{Int}   # lookup the row block given the row
    colblocks::Vector{Int}   # lookup col row block given the row

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



## View

typealias SubBandedBlockSubBlock{T,BBM<:BandedBlockMatrix,II<:Union{Block,SubBlock},JJ<:Union{Block,SubBlock}} SubArray{T,2,BBM,Tuple{II,JJ},false}
typealias SubBandedBlockRange{T,BBM<:BandedBlockMatrix} SubArray{T,2,BBM,Tuple{StepRange{Block,Block},StepRange{Block,Block}},false}
typealias StridedMatrix2{T,A<:Union{DenseArray,Base.StridedReshapedArray,BandedBlockMatrix},I<:Tuple{Vararg{Union{Base.RangeIndex, Base.AbstractCartesianIndex, Block,SubBlock}}}}  Union{DenseArray{T,2}, SubArray{T,2,A,I}, Base.StridedReshapedArray{T,2}}


subblocksize(A::BandedBlockMatrix,K::Block,J::Block) =
    (A.rows[K.K],A.cols[J.K])
subblocksize(A::BandedBlockMatrix,K::Block,J::SubBlock) =
    (A.rows[K.K],length(J.sub))
subblocksize(A::BandedBlockMatrix,K::SubBlock,J::Block) =
    (length(K.sub),A.cols[J.K])
subblocksize(A::BandedBlockMatrix,K::SubBlock,J::SubBlock) =
    (length(K.sub),length(J.sub))

subblocksize(A::BandedBlockMatrix,K::StepRange{Block,Block},J::StepRange{Block,Block}) =
    (sum(A.rows[Int.(K)]),sum(A.rows[Int.(J)]))

Base.view(A::BandedBlockMatrix,K::Union{Block,SubBlock,StepRange{Block,Block}},J::Union{Block,SubBlock,StepRange{Block,Block}}) =
    SubArray(A, (K,J), subblocksize(A,K,J))

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

rowblocklengths(A::SubBandedBlockRange) = parent(A).rows[parentindexes(A)[1]]
colblocklengths(A::SubBandedBlockRange) = parent(A).cols[parentindexes(A)[2]]

Base.view(A::SubBandedBlockRange,K::Union{Block,SubBlock,StepRange{Block,Block}},J::Union{Block,SubBlock,StepRange{Block,Block}}) =
    view(parent(A),parentindexes(A)[1][Int.(K)],parentindexes(A)[2][Int.(J)])


function Base.indices(S::Union{SubBandedBlockSubBlock,SubBandedBlockRange})
    sz = subblocksize(parent(S),parentindexes(S)...)
    (Base.OneTo(sz[1]),
     Base.OneTo(sz[2]))
 end

parentblock(S::SubBandedBlockSubBlock) =
     view(parent(S),block(parentindexes(S)[1]),block(parentindexes(S)[2]))


# returns a view of the data corresponding to the block
function blockview{T,U,V}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{Block,Block},false})
    A = parent(S)
    K,J = parentindexes(S)
    st = A.blockstart[block(K).K,block(J).K]
    sz = size(S)
    reshape(view(A.data,st+1:st+sz[1]*sz[2]),sz[1],sz[2])
end

# returns a view of the data
dataview{T,U,V}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{Block,Block},false}) =
    blockview(S)
dataview{T,U,V,II}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{SubBlock{II},Block},false}) =
    view(blockview(parentblock(S)),parentindexes(S)[1].sub,:)
dataview{T,U,V,JJ}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{Block,SubBlock{JJ}},false}) =
    view(blockview(parentblock(S)),:,parentindexes(S)[2].sub)
dataview{T,U,V,II,JJ}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{SubBlock{II},SubBlock{JJ}},false}) =
    view(blockview(parentblock(S)),parentindexes(S)[1].sub,parentindexes(S)[2].sub)

getindex{N}(S::SubBandedBlockSubBlock, I::Vararg{Real,N}) =
    dataview(S)[I...]

function getindex(S::SubBandedBlockRange, k::Integer, j::Integer)
    KR,JR = parentindexes(S)
    @assert step(KR) == step(JR) == Block(1)
    A = parent(S)
    A[k + sum(A.rows[1:first(KR).K-1]),j + sum(A.rows[1:first(JR).K-1])]
end

function setindex!{N}(S::SubBandedBlockSubBlock, v, I::Vararg{Real,N})
    dataview(S)[I...] = v
end

function setindex!(S::SubBandedBlockRange, v, k::Integer, j::Integer)
    KR,JR = parentindexes(S)
    @assert step(KR) == step(JR) == Block(1)
    A = parent(S)
    A[k + sum(A.rows[1:first(KR).K-1]),j + sum(A.rows[1:first(JR).K-1])] = v
end


Base.strides(S::SubBandedBlockSubBlock) =
    (1,parent(S).rows[block(parentindexes(S)[1]).K])

function Base.pointer{T<:BlasFloat,U,V}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{Block,Block},false})
    A = parent(S)
    K,J = parentindexes(S)
    pointer(A.data)+A.blockstart[K.K,J.K]*sizeof(T)
end


Base.pointer{T<:BlasFloat,U,V,JJ}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{Block,SubBlock{JJ}},false}) =
    pointer(parentblock(S))

Base.pointer{T<:BlasFloat,U,V,II}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{SubBlock{II},Block},false}) =
    pointer(parentblock(S)) + sizeof(T)*(first(parentindexes(S)[1].sub)-1)

Base.pointer{T<:BlasFloat,U,V,II,JJ}(S::SubArray{T,2,BandedBlockMatrix{T,U,V},Tuple{SubBlock{II},SubBlock{JJ}},false}) =
    pointer(parentblock(S)) + sizeof(T)*(first(parentindexes(S)[1].sub)-1) +
        sizeof(T)*stride(S,2)*(first(parentindexes(S)[2].sub)-1)


## algebra

αA_mul_B_plus_βC!(α,A::SubBandedBlockSubBlock,x::AbstractVector,β,y::AbstractVector) =
    gemv!('N',α,A,x,β,y)
αA_mul_B_plus_βC!(α,A::StridedMatrix2,B::StridedMatrix2,β,C::StridedMatrix2) =
    gemm!('N','N',α,A,B,β,C)

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
        BandedBlockMatrix(copy(vec(B.data)),0,0,ones(Int,size(B,1)),ones(Int,size(B,2)))
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

function bandedblock_rectblocks_trtrs!{T}(R::BandedBlockMatrix{T},b::Vector)
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


function trtrs!{T}(A::BandedBlockMatrix{T},u::Matrix)
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
