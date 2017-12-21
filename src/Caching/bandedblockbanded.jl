
function CachedOperator(::Type{BandedBlockBandedMatrix}, op::Operator)
    l,u = blockbandwidths(op)
    λ,μ = subblockbandwidths(op)
    data = BandedBlockBandedMatrix{eltype(op)}(uninitialized,
        (blocklengths(rangespace(op))[1:0],blocklengths(domainspace(op))[1:0]),
        (l,u), (λ,μ))

    CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),false)
end

# Grow cached operator
#
function resizedata!(B::CachedOperator{T,BandedBlockBandedMatrix{T}}, ::Colon, col::Integer) where {T<:Number}
    if col > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if col > B.datasize[2]
        l,u,λ,μ = B.data.l,B.data.u,B.data.λ,B.data.μ
        J = Int(block(domainspace(B),col))

        rows = blocklengths(rangespace(B.op))[1:J+l]
        cols = blocklengths(domainspace(B.op))[1:J]

        bs = BandedBlockBandedSizes(rows, cols, l, u, λ, μ)

        # B.data = _BandedBlockBandedMatrix(PseudoBlockArray

        blocks = PseudoBlockArray(pad(B.data.data.blocks,:,(l+u+1)*sum(cols)), bs.data_block_sizes)
        B.data = _BandedBlockBandedMatrix(blocks, bs)

        jr=B.datasize[2]+1:col
        kr=colstart(B.data,jr[1]):colstop(B.data,jr[end])

        isempty(kr) || BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (last(kr),col)
    end

    B
end

function resizedata!(B::CachedOperator{T,BandedBlockBandedMatrix{T}},n::Integer,m::Integer) where {T<:Number}
    resizedata!(B,:,m)
    if n < B.datasize[1]
        return B
    end
    
    l,u,λ,μ = B.data.l,B.data.u,B.data.λ,B.data.μ

    # make sure we have enough rows
    K = Int(block(rangespace(B),n))
    rows = blocklengths(rangespace(B.op))[1:K]

    b_bs = BlockSizes(BlockArrays._cumul_vec(rows), B.data.block_sizes.cumul_sizes[2])

    bs = BandedBlockBandedSizes(b_bs, B.data.block_sizes.data_block_sizes,
                                l, u, λ, μ)

    B.data.rows=rows
    B.data.rowblocks=blocklookup(rows)
    B.datasize=(n,B.datasize[2])

    B
end
