function CachedOperator(::Type{BandedBlockBandedMatrix},op::Operator)
    l,u=blockbandwidths(op)
    λ,μ=subblockbandwidths(op)
    data=BandedBlockBandedMatrix(eltype(op),l,u,λ,μ,
                            ApproxFun.blocklengths(rangespace(op))[1:0],
                            ApproxFun.blocklengths(domainspace(op))[1:0])
    CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),false)
end


# Grow cached operator
#
function resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockBandedMatrix{T,RI,DI}},::Colon,col::Integer)
    if col > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if col > B.datasize[2]
        l=B.data.l; u=B.data.u
        J=block(domainspace(B),col)

        rows=blocklengths(rangespace(B.op))[1:J+l]
        cols=blocklengths(domainspace(B.op))[1:J]

        B.data.data=pad(B.data.data,:,(l+u+1)*sum(cols))
        B.data.rows=rows
        B.data.rowblocks=blocklookup(rows)
        B.data.cols=cols
        B.data.colblocks=blocklookup(cols)

        jr=B.datasize[2]+1:col
        kr=colstart(B.data,jr[1]):colstop(B.data,jr[end])
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (kr[end],col)
    end

    B
end

resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockBandedMatrix{T,RI,DI}},n::Integer,m::Integer) =
    resizedata!(B,:,m)
