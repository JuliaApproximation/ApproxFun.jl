


# diagblockshift gives the shift for the diagonal block of an operator


diagblockshift(a,b) = error("Developer: Not implemented for blocklengths $a, $b")

diagblockshift(a::UnitCount,b::UnitCount) = b.start-a.start

diagblockshift(op::Operator) =
    diagblockshift(blocklengths(domainspace(op)),blocklengths(rangespace(op)))


function CachedOperator(::Type{BandedBlockMatrix},op::Operator;padding::Bool=false)
    l,u=blockbandwidths(op)
    padding && (u+=l+diagblockshift(op))
    data=BandedBlockMatrix(eltype(op),l,u,1:0,1:0)  # TODO: type of rows/cols
    CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),padding)
end



## colstop, etc.

for OP in (:colstop,:colstart,:rowstop,:rowstart)
    @eval begin
        $OP{T,RI,DI,M,DS,RS}(D::CachedOperator{T,BandedBlockMatrix{T,RI,DI},M,DS,RS,Tuple{Infinity{Bool},Infinity{Bool}}},
                             k::Integer) = $OP(D.op,k)
        $OP{T,RI,DI,M,DS,RS}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI},M,DS,RS},k::Integer) = $OP(B.op,k)
    end
end



## Grow cached operator
#
function resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},::Colon,col::Integer)
    if col > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if col > B.datasize[2]
        l=B.data.l; u=B.data.u
        J=block(domainspace(B),col)

        rows=blocklengths(rangespace(B.op))[1:J+l]
        cols=blocklengths(domainspace(B.op))[1:J]

        pad!(B.data.data,bbm_numentries(rows,cols,l,u))
        B.data.rows=rows
        B.data.rowblocks=blocklookup(rows)
        B.data.cols=cols
        B.data.colblocks=blocklookup(cols)
        B.data.blockstart=bbm_blockstarts(rows,cols,l,u)

        jr=B.datasize[2]+1:col
        kr=colstart(B.data,jr[1]):colstop(B.data,jr[end])
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (kr[end],col)
    end

    B
end

resizedata!{T<:Number,RI,DI}(B::CachedOperator{T,BandedBlockMatrix{T,RI,DI}},n::Integer,m::Integer) =
    resizedata!(B,:,m)


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



function resizedata!{T<:BlasFloat,MM,DS,RS,DDS,RRS,BI}(QR::QROperator{CachedOperator{T,BandedBlockMatrix{T,DDS,RRS},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data
    # last block, convoluted def to match blockbandedmatrix
    J_col = block(domainspace(MO),col)
    K_end = blockcolstop(MO,J_col)  # last row block in last column
    J_end = blockrowstop(MO,K_end)  # QR will affect up to this column
    rs=blockstop(rangespace(MO),J_end)  # we need to resize up this column
    sz=sizeof(T)

    if rs ≥ MO.datasize[2]
        # add columns up to column rs, which is last column affected by QR
        resizedata!(MO,:,rs)
        R=MO.data
    end

    if col > size(W,2)
        # resize Householder matrix
        m=size(W,2)
        resize!(W.cols,2col+1)

        for j=m+1:2col
            cs = blockstop(rangespace(MO),blockcolstop(MO,block(domainspace(MO),j)))
            W.cols[j+1]=W.cols[j] + cs-j+1
            W.m=max(W.m,cs-j+1)
        end

        resize!(W.data,W.cols[end]-1)
    end

    w=pointer(W.data)
    r=pointer(R.data)

    for k=QR.ncols+1:col
        J1=R.colblocks[k]
        CS=blockcolstop(R,J1)

        wp=w+sz*(W.cols[k]-1)          # k-th column of W



        K1=R.rowblocks[k]  # diagonal block
        br=blockrows(R,K1)  # the rows of the diagonal block
        bc=blockcols(R,J1) # the cols of the diagonal block

        # copy into W
        M1 = br[end]-k+1 # the number of entries in first block
        kshft = k-br[1]
        nrows1 = R.rows[K1]
        jshft = (k-bc[1])*nrows1  # shift by each column we are to the right
        BLAS.blascopy!(M1,r+sz*(R.blockstart[K1,J1]+kshft + jshft),1,wp,1)

        M = M1 # number of entries so far copied


        for K=K1+1:CS
            jshft = (k-bc[1])*R.rows[K]  # shift by each column we are to the right
            BLAS.blascopy!(R.rows[K],r+sz*(R.blockstart[K,J1] + jshft),1,wp+sz*M,1)
            M += R.rows[K]    # copy all rows in K-th block
        end

        # M is now total entries in W
        W.data[W.cols[k]] += flipsign(BLAS.nrm2(M,wp,1),W.data[W.cols[k]])
        normalize!(M,wp)

        # first block
        # scale banded entries
        BRS1=blockrowstop(R,K1)
        for J=J1:BRS1
            for j=(J==J1?k-bc[1]+1:1):R.cols[J]  # only do partial columns for first block
                jshft = (j-1)*nrows1
                dt=dot(M1,wp,1,r+sz*(R.blockstart[K1,J]+kshft +jshft),1)
                M=M1
                for K=K1+1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    dt+=dot(R.rows[K],wp+sz*M,1,r+sz*(R.blockstart[K,J] +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end

                jshft = (j-1)*nrows1
                BLAS.axpy!(M1,-2*dt,wp,1,r+sz*(R.blockstart[K1,J]+kshft +jshft),1)
                M=M1
                for K=K1+1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    BLAS.axpy!(R.rows[K],-2*dt,wp+sz*M,1,r+sz*(R.blockstart[K,J] +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
            end
        end


        # now do the blocks where we have zeros
        for J=BRS1+1:blockrowstop(R,CS)
            for j=1:R.cols[J]  # only do partial columns for first block
                dt=zero(T)
                Mpre=M1 + sum(R.rows[K1+1:K1+J-BRS1-1]) # number of rows in zero blocks
                M=Mpre
                for K=K1+J-BRS1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    dt+=dot(R.rows[K],wp+sz*M,1,r+sz*(R.blockstart[K,J] +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
                M=Mpre
                for K=K1+J-BRS1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    BLAS.axpy!(R.rows[K],-2*dt,wp+sz*M,1,r+sz*(R.blockstart[K,J] +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
            end
        end
    end
    QR.ncols=col
    QR
end
