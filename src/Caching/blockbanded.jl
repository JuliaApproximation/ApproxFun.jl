


# # default copy is to loop through
# # override this for most operators.
function default_BlockBandedMatrix(S::Operator)
    ret = BlockBandedMatrix(Zeros, S)

    @inbounds for J=Block(1):Block(nblocks(ret,2)), K=blockcolrange(ret,Int(J))
        ret[K,J] = view(S,K,J)
    end
    ret
end
#


# diagblockshift gives the shift for the diagonal block of an operator
# that is, trace an operator down the diagonal.  What blocks correspond to the
# block diagonal?
# this is used to determine how many blocks to pad in the QR decomposition, as
# every lower block gets added to the upper daigonal

diagblockshift(a,b) = error("Developer: Not implemented for blocklengths $a, $b")

function diagblockshift(a::AbstractRange, b::AbstractRange)
    @assert step(a) == step(b)
    b.start-a.start
end
diagblockshift(a::AbstractUnitRange, b::AbstractUnitRange) = b.start-a.start

#TODO: generalize
function diagblockshift(a::AbstractFill{Int},b::AbstractFill{Int})
    @assert getindex_value(a) == getindex_value(b)
    0
end

diagblockshift(a::AbstractFill{Int},b::Vcat{Int,1,<:Tuple{V1,<:AbstractFill{Int}}}) where {V1 <: AbstractVector{Int}} =
    max(0,-diagblockshift(b,a))


function diagblockshift(a::Vcat{Int,1,<:Tuple{V1,<:AbstractFill{Int}}},b::AbstractFill{Int}) where V1 <: AbstractVector{Int}
    @assert getindex_value(a.it[end]) == getindex_value(b)
    isempty(a.it[1]) && return diagblockshift(a.it[2],b)
    a1, b1 = a[1],b[1]
    a1 == b1 && return diagblockshift(Vcat(a.it[1][2:end],a.it[2]),b)
    a1 >  b1 && length(a.it[1]) == 1 && return 0
    a1 >  b1 && return max(0,-1+diagblockshift(flatten(([a1-b1;a.it[1][2:end]],a.it[2]),b)))
    a1 <  b1 && length(a.it[1]) == 1 && return 1
    # a1 <  b1 &&
    return 1+diagblockshift(Vcat(a.it[1][2:end],a.it[2]),Vcat([b1-a1],b))
end

function diagblockshift(a::Vcat{Int,1,<:Tuple{V1,<:AbstractFill{Int}}},
                        b::Vcat{Int,1,<:Tuple{V2,<:AbstractFill{Int}}}) where {V1 <: AbstractVector{Int},V2 <: AbstractVector{Int}}
    isempty(a.it[1]) && return diagblockshift(a.it[2],b)
    isempty(b.it[1]) && return diagblockshift(a,b.it[2])
    a1, b1 = a[1],b[1]
    a1 == b1 && return diagblockshift(Vcat(a.it[1][2:end],a.it[2]),Vcat(b.it[1][2:end],b.it[2]))
    a1 >  b1 && return max(0,-1+diagblockshift(Vcat([a1-b1;a.it[1][2:end]],a.it[2]),
                                               Vcat(b.it[1][2:end],b.it[2])))
    # a1 <  b1 &&
    return 1+diagblockshift(Vcat(a.it[1][2:end],a.it[2]),Vcat([b1-a1;b.it[1][2:end]],b.it[2]))
end


diagblockshift(op::Operator) =
    diagblockshift(blocklengths(domainspace(op)),blocklengths(rangespace(op)))


function CachedOperator(::Type{BlockBandedMatrix},op::Operator;padding::Bool=false)
    l,u=blockbandwidths(op)
    padding && (u+=l+diagblockshift(op))
    data=BlockBandedMatrix{eltype(op)}(undef,
                                        (Vector{Int}(), Vector{Int}()),
                                        (l,u))
    CachedOperator(op,data,(0,0),domainspace(op),rangespace(op),(-l,u),padding)
end





## Grow cached operator
#
function resizedata!(B::CachedOperator{T,BlockBandedMatrix{T}},::Colon,col::Integer) where {T<:Number}
    if col > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if col > B.datasize[2]
        if B.datasize[2] == 0
            datablocksize = Block(0)
        else
            datablocksize = block(domainspace(B),B.datasize[2])
            bs = blockstop(domainspace(B),datablocksize)
            if bs ≠ B.datasize[2]
                error("Internal Error: $(B.datasize) is not lined up with the block $datablocksize as the last column doesn't end at $bs")
            end
        end


        l,u = blockbandwidths(B.data)
        J = block(domainspace(B),col)
        col = blockstop(domainspace(B),J)  # pad to end of block

        rows = blocklengths(rangespace(B.op))[1:Int(J)+l]
        cols = blocklengths(domainspace(B.op))[1:Int(J)]

        b_size = BlockBandedSizes(Vector{Int}(rows), Vector{Int}(cols), l, u)

        resize!(B.data.data, bb_numentries(b_size))
        B.data = _BlockBandedMatrix(B.data.data, b_size)

        JR = datablocksize+1:J
        KR=blockcolstart(B.data,first(JR)):blockcolstop(B.data,last(JR))
        copyto!(view(B.data,KR,JR), view(B.op,KR,JR))

        B.datasize = (blockstop(rangespace(B),last(KR)),col)
    end

    B
end

function resizedata!(B::CachedOperator{T,BlockBandedMatrix{T}},n::Integer,m::Integer) where {T<:Number}
    N = block(rangespace(B),n)
    m̃ = blockstart(domainspace(B),N)
    resizedata!(B,:,max(m,m̃))
end


## QR
# we use a RaggedMatrix to represent the growing lengths of the
# householder reflections
QROperator(R::CachedOperator{T,BlockBandedMatrix{T}}) where {T} =
    QROperator(R,RaggedMatrix{T}(undef, 0, Int[]),0)


# function resizedata!(QR::QROperator{CachedOperator{T,BlockBandedMatrix{T},
#                                           MM,DS,RS,BI}},
#                      ::Colon, col) where {T,MM,DS,RS,BI}
#     if col ≤ QR.ncols
#         return QR
#     end
#
#     MO=QR.R_cache
#     W=QR.H
#
#     R=MO.data
#     cs=colstop(MO,col)
#
#     if cs ≥ MO.datasize[1]
#         resizedata!(MO,cs+100,:)  # add 100 rows
#         R=MO.data
#     end
#
#     if col > size(W,2)
#         m=size(W,2)
#         resize!(W.cols,2col+1)
#
#         for j=m+1:2col
#             cs=colstop(MO,j)
#             W.cols[j+1]=W.cols[j] + cs-j+1
#             W.m=max(W.m,cs-j+1)
#         end
#
#         resize!(W.data,W.cols[end]-1)
#     end
#
#     for k=QR.ncols+1:col
#         cs = colstop(R,k)
#         W[1:cs-k+1,k] = view(R,k:cs,k) # diagonal and below
#         wp=view(W,1:cs-k+1,k)
#         W[1,k]+= flipsign(norm(wp),W[1,k])
#         normalize!(wp)
#
#         # scale banded entries
#         for j=k:rowstop(R,k)
#             v=view(R,k:cs,j)
#             dt=dot(wp,v)
#             LinearAlgebra.axpy!(-2*dt,wp,v)
#         end
#
#         # scale banded/filled entries
#         for j=rowstop(R,k)+1:rowstop(R,cs)
#             csrt=colstart(R,j)
#             v=view(R,csrt:cs,j)  # shift down each time
#             wp2=view(wp,csrt-k+1:cs-k+1)
#             dt=dot(wp2,v)
#             LinearAlgebra.axpy!(-2*dt,wp2,v)
#         end
#     end
#     QR.ncols=col
#     QR
# end

# always resize by column
resizedata!(QR::QROperator{CachedOperator{T,BlockBandedMatrix{T},
                               MM,DS,RS,BI}},
                     ::Colon, col::Int) where {T,MM,DS,RS,BI} =
    resizedata!(QR, :, block(domainspace(QR.R_cache),col))

function resizedata!(QR::QROperator{CachedOperator{T,BlockBandedMatrix{T},
                               MM,DS,RS,BI}},
                     ::Colon, COL::Block) where {T<:BlasFloat,MM,DS,RS,BI}
     MO = QR.R_cache
     W = QR.H
     R = MO.data

     l,u = R.l,R.u

     ds = domainspace(QR)
     col = blockstop(ds, COL)  # last column

     if col ≤ QR.ncols
         return QR
     end

     J_start = Int(block(ds, QR.ncols+1)) # first block-column not factorized
     @assert blockstart(ds, J_start) == QR.ncols+1  # we want to make sure we haven't partially factorized a column

     # last block, convoluted def to match blockbandedmatrix
     J_col = Int(COL)
     K_end = Int(blockcolstop(MO, J_col))  # last row block in last column
     J_end = Int(blockrowstop(MO, K_end))  # QR will affect up to this column
     j_end = blockstop(domainspace(MO), J_end)  # we need to resize up this column
     sz = sizeof(T)

     if j_end ≥ MO.datasize[2]
         # add columns up to column rs, which is last column affected by QR
         resizedata!(MO, :, j_end)
         R = MO.data
     end

     if col > size(W,2)
         # resize Householder matrix
         m = size(W,2)
         resize!(W.cols, 2col+1)

         for j=m+1:2col
             cs = blockstop(rangespace(MO), blockcolstop(MO, block(domainspace(MO),j)))
             W.cols[j+1]=W.cols[j] + cs-j+1
             W.m=max(W.m,cs-j+1)
         end

         resize!(W.data, W.cols[end]-1)
     end

     w = pointer(W.data)
     r = pointer(R.data)

     bs = R.block_sizes

     for j =QR.ncols+1 : col   # first column of block
         bi = global2blockindex(bs, (j, j)) # converts from global indices to block indices
         K1, J1 = bi.I  # this is the diagonal block corresponding to j
         κ, ξ = bi.α

         st = bs.block_strides[J1]  # the stride of the matrix
         shft = bs.block_starts[K1,J1]-1 + st*(ξ-1) + κ-1 # the index of the pointer to the j, j entry


         K_CS = Int(blockcolstop(R,J1)) # last row in J-th blockcolumn
         k_end = globalrange(bs, (K_CS, J1))[1][end]


         w_j = W.cols[j]  # the data index  for the j-th column of W
         wp = w+sz*(w_j-1)          # j-th column of W

         M = k_end - j + 1 # the number of entries we are diagonalizing. we know the stride tells us the total number of rows

         BLAS.blascopy!(M, r+sz*shft, 1, wp, 1) # copy the column into W


         # we need to scale the first entry and then normalize
         W.data[w_j] += flipsign(BLAS.nrm2(M,wp,1), W.data[w_j])
         normalize!(M, wp)

         # scale rest of columns in first block
         # for ξ_2 = 2:

         for ξ_2 = ξ:blocksize(bs, 2, J1)
             # we now apply I-2v*v' in place
             r_sh = r+sz*(shft + st*(ξ_2-ξ)) # the pointer the (j,ξ_2)-th entry
             dt = BandedMatrices.dot(M, wp, 1, r_sh, 1)
             BLAS.axpy!(M, -2*dt, wp, 1, r_sh ,1)
         end

         for J = J1+1:min(K1+u,J_end)
             st = bs.block_strides[J]
             shft = bs.block_starts[K1,J] + κ-2 # the index of the pointer to the j, j entry
             for ξ_2 = 1:blocksize(bs.block_sizes, 2, J)
                 # we now apply I-2v*v' in place
                 # r_sh = r+sz*(shft + st*(ξ_2-1)) # the pointer the (j,ξ_2)-th entry

                 # TODO: remove these debugging statement
                 # @assert w_j-1 + M ≤ length(W.data)
                 # @assert shft + st*(ξ_2-1) + M ≤ length(R.data)
                 # @assert 0 ≤ w_j-1
                 # if ! (0 ≤ shft + st*(ξ_2-1))
                 #     @show shft, st, ξ_2, l, u
                 #     @show κ, bs.block_starts[K1,J]
                 #     @show K1, J
                 #     @show MO.op
                 # end
                 # dt = dot(M, wp, 1, r_sh, 1)
                 # BLAS.axpy!(M, -2*dt, wp, 1, r_sh ,1)

                 dt = dot(view(W.data, w_j:w_j+M-1) ,
                                    view(R.data, shft + st*(ξ_2-1) +1:shft + st*(ξ_2-1) +M))
                 BLAS.axpy!(-2*dt, view(W.data, w_j:w_j+M-1) ,
                                    view(R.data, shft + st*(ξ_2-1) +1:shft + st*(ξ_2-1) +M))
             end
         end
     end

    QR.ncols=col
    QR
end
