
function CachedOperator(::Type{BandedBlockMatrix},op::Operator;padding::Bool=false)
    l,u=blockbandwidths(op)
    padding=false
    padding && (u+=l)
    data=BandedBlockMatrix(eltype(op),l,u,1:0,1:0)  # TODO: type of rows/cols
    CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),padding)
end



## colstop, etc.

for OP in (:colstop,:colstart,:rowstop,:rowstart)
    @eval $OP{T,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},k::Integer) = $OP(B.op,k)
end



## Grow cached operator
#
function resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},n::Integer,::Colon)
    if n > size(B,1)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if n > B.datasize[1]
        l=B.data.l; u=B.data.u
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


## QR
# we use a RaggedMatrix to represent the growing lengths of the
# householder reflections
QROperator{T,DDS,RRS}(R::CachedOperator{T,BandedBlockMatrix{T,DDS,RRS}}) =
    QROperator(R,RaggedMatrix(T,0,Int[]),0)


function resizedata!{T,MM,DS,RS,DDS,RRS,BI}(QR::QROperator{CachedOperator{T,BandedBlockMatrix{T,DDS,RRS},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data
    cs=colstop(MO,col)

    if cs ≥ MO.datasize[1]
        resizedata!(MO,cs+100,:)  # add 100 rows
        R=MO.data
    end

    if col > size(W,2)
        m=size(W,2)
        resize!(W.cols,2col+1)

        for j=m+1:2col
            cs=colstop(MO,j)
            W.cols[j+1]=W.cols[j] + cs-j+1
            W.m=max(W.m,cs-j+1)
        end

        resize!(W.data,W.cols[end]-1)
    end

    for k=QR.ncols+1:col
        cs = colstop(R,k)
        W[1:cs-k+1,k] = view(R,k:cs,k) # diagonal and below
        wp=view(W,1:cs-k+1,k)
        W[1,k]+= flipsign(norm(wp),W[1,k])
        normalize!(wp)

        # scale banded entries
        for j=k:rowstop(R,k)
            v=view(R,k:cs,j)
            dt=dot(wp,v)
            Base.axpy!(-2*dt,wp,v)
        end

        # scale banded/filled entries
        for j=rowstop(R,k)+1:rowstop(R,cs)
            csrt=colstart(R,j)
            v=view(R,csrt:cs,j)  # shift down each time
            wp2=view(wp,csrt-k+1:cs-k+1)
            dt=dot(wp2,v)
            Base.axpy!(-2*dt,wp2,v)
        end
    end
    QR.ncols=col
    QR
end
