
## Grow cached operator
#
function resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},n::Integer,::Colon)
    if n > size(B,1)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if n > B.datasize[1]
        K=block(rangespace(B),n)

        rows=blocklengths(rangespace(B.op))[1:K]
        cols=blocklengths(domainspace(B.op))[1:K+B.data.u]

        pad!(B.data.data,bbm_numentries(rows,cols,l,u))
        B.data.rows=rows
        B.data.rowblocks=blocklookup(rows)
        B.data.cols=cols
        B.data.colblocks=blocklookup(cols)
        B.data.blockstart=bbm_blockstarts(rows,cols,l,u)

        kr=B.datasize[1]+1:n
        jr=max(blockstart(domainspace(B),block(rangespace(B),B.datasize[1]+1)-B.data.l),1):min(blockstop(domainspace(B),K+B.data.u),size(B,2))
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (n,jr[end])
    end

    B
end

resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},n::Integer,m::Integer) =
    resizedata!(B,n,:)
