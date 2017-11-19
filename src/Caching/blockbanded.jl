


# # default copy is to loop through
# # override this for most operators.
function default_blockbandedmatrix(S::Operator)
    ret=BlockBandedMatrix(eltype(S),blockbandwidth(S,1),blockbandwidth(S,2),
            blocklengths(rangespace(S)),blocklengths(domainspace(S)))

    @inbounds for J=Block(1):Block(blocksize(ret,2)),K=blockcolrange(ret,J)
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

function diagblockshift(a::AbstractCount,b::AbstractCount)
    @assert step(a) == step(b)
    b.start-a.start
end
diagblockshift(a::UnitCount,b::UnitCount) = b.start-a.start

#TODO: generalize
function diagblockshift(a::Repeated{Int},b::Repeated{Int})
    @assert a.x == b.x
    0
end

diagblockshift(a::Repeated{Int},b::Flatten{Tuple{V1,Repeated{Int}}}) where {V1 <: AbstractVector{Int}} =
    max(0,-diagblockshift(b,a))


function diagblockshift(a::Flatten{Tuple{V1,Repeated{Int}}},b::Repeated{Int}) where V1 <: AbstractVector{Int}
    @assert a.it[end].x == b.x
    isempty(a.it[1]) && return diagblockshift(a.it[2],b)
    a1, b1 = a[1],b[1]
    a1 == b1 && return diagblockshift(flatten((a.it[1][2:end],a.it[2])),b)
    a1 >  b1 && length(a.it[1]) == 1 && return 0
    a1 >  b1 && return max(0,-1+diagblockshift(flatten(([a1-b1;a.it[1][2:end]],a.it[2]),b)))
    a1 <  b1 && length(a.it[1]) == 1 && return 1
    # a1 <  b1 &&
    return 1+diagblockshift(flatten((a.it[1][2:end],a.it[2])),flatten(([b1-a1],b)))
end

function diagblockshift(a::Flatten{Tuple{V1,Repeated{Int}}},
                        b::Flatten{Tuple{V2,Repeated{Int}}}) where {V1 <: AbstractVector{Int},V2 <: AbstractVector{Int}}
    isempty(a.it[1]) && return diagblockshift(a.it[2],b)
    isempty(b.it[1]) && return diagblockshift(a,b.it[2])
    a1, b1 = a[1],b[1]
    a1 == b1 && return diagblockshift(flatten((a.it[1][2:end],a.it[2])),flatten((b.it[1][2:end],b.it[2])))
    a1 >  b1 && return max(0,-1+diagblockshift(flatten(([a1-b1;a.it[1][2:end]],a.it[2])),
                                               flatten((b.it[1][2:end],b.it[2]))))
    # a1 <  b1 &&
    return 1+diagblockshift(flatten((a.it[1][2:end],a.it[2])),flatten(([b1-a1;b.it[1][2:end]],b.it[2])))
end


diagblockshift(op::Operator) =
    diagblockshift(blocklengths(domainspace(op)),blocklengths(rangespace(op)))


function CachedOperator(::Type{BlockBandedMatrix},op::Operator;padding::Bool=false)
    l,u=blockbandwidths(op)
    padding && (u+=l+diagblockshift(op))
    data=BlockBandedMatrix(eltype(op),l,u,
                            blocklengths(rangespace(op))[1:0],
                            blocklengths(domainspace(op))[1:0])
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
            datablocksize = 0
        else
            datablocksize = block(domainspace(B),B.datasize[2])
            bs =  blockstop(domainspace(B),datablocksize)
            if bs ≠ B.datasize[2]
                error("Developer: $(B.datasize) is not lined up with the block $datablocksize as the last column doesn't end at $bs")
            end
        end


        l=B.data.l; u=B.data.u
        J=block(domainspace(B),col)
        col = blockstop(domainspace(B),J)  # pad to end of block

        rows = blocklengths(rangespace(B.op))[1:J.n[1]+l]
        cols = blocklengths(domainspace(B.op))[1:J.n[1]]

        pad!(B.data.data,bbm_numentries(rows,cols,l,u))
        B.data.rows = rows
        B.data.rowblocks = blocklookup(rows)
        B.data.cols = cols
        B.data.colblocks = blocklookup(cols)
        B.data.blockstart = bbm_blockstarts(rows,cols,l,u)

        JR=datablocksize+1:J
        KR=blockcolstart(B.data,JR[1]):blockcolstop(B.data,JR[end])
        BLAS.axpy!(1.0,view(B.op,KR,JR),view(B.data,KR,JR))

        B.datasize = (blockstop(rangespace(B),KR[end]),col)
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
    QROperator(R,RaggedMatrix(T,0,Int[]),0)


function resizedata!(QR::QROperator{CachedOperator{T,BlockBandedMatrix{T},
                                          MM,DS,RS,BI}},
 ::Colon,col) where {T,MM,DS,RS,BI}
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



function resizedata!(QR::QROperator{CachedOperator{T,BlockBandedMatrix{T},
                               MM,DS,RS,BI}},
::Colon,col) where {T<:BlasFloat,MM,DS,RS,BI}
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data
    # last block, convoluted def to match blockbandedmatrix
    J_col = block(domainspace(MO),col).n[1]
    K_end = blockcolstop(MO,J_col).n[1]  # last row block in last column
    J_end = blockrowstop(MO,K_end).n[1]  # QR will affect up to this column
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

    @inbounds for k=QR.ncols+1:col
        J1=R.colblocks[k]
        CS=blockcolstop(R,J1).n[1]

        wp=w+sz*(W.cols[k]-1)          # k-th column of W



        K1=R.rowblocks[k]  # diagonal block
        br=blockrows(R,K1)  # the rows of the diagonal block
        bc=blockcols(R,J1) # the cols of the diagonal block

        # copy into W
        M1 = br[end]-k+1 # the number of entries in first block
        kshft = k-br[1]
        nrows1 = R.rows[K1]
        jshft = (k-bc[1])*nrows1  # shift by each column we are to the right
        BLAS.blascopy!(M1,r+sz*(inbands_getindex(R.blockstart,K1,J1)+kshft + jshft),1,wp,1)

        M = M1 # number of entries so far copied


        for K=K1+1:CS
            jshft = (k-bc[1])*R.rows[K]  # shift by each column we are to the right
            BLAS.blascopy!(R.rows[K],r+sz*(inbands_getindex(R.blockstart,K,J1) + jshft),1,wp+sz*M,1)
            M += R.rows[K]    # copy all rows in K-th block
        end

        # M is now total entries in W
        W.data[W.cols[k]] += flipsign(BLAS.nrm2(M,wp,1),W.data[W.cols[k]])
        normalize!(M,wp)

        # first block
        # scale banded entries
        BRS1=blockrowstop(R,K1).n[1]
        @inbounds for J=J1:BRS1
            for j = (J==J1 ? k-bc[1]+1 : 1):R.cols[J]  # only do partial columns for first block
                jshft = (j-1)*nrows1
                dt=dot(M1,wp,1,r+sz*(inbands_getindex(R.blockstart,K1,J)+kshft +jshft),1)
                M=M1
                for K=K1+1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    dt+=dot(R.rows[K],wp+sz*M,1,r+sz*(inbands_getindex(R.blockstart,K,J) +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end

                jshft = (j-1)*nrows1
                BLAS.axpy!(M1,-2*dt,wp,1,r+sz*(inbands_getindex(R.blockstart,K1,J)+kshft +jshft),1)
                M=M1
                for K=K1+1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    BLAS.axpy!(R.rows[K],-2*dt,wp+sz*M,1,r+sz*(inbands_getindex(R.blockstart,K,J) +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
            end
        end


        # now do the blocks where we have zeros
        @inbounds for J=BRS1+1:blockrowstop(R,CS).n[1]
            for j=1:R.cols[J]  # only do partial columns for first block
                dt=zero(T)
                Mpre=M1 + sum(R.rows[K1+1:K1+J-BRS1-1]) # number of rows in zero blocks
                M=Mpre
                for K=K1+J-BRS1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    dt+=dot(R.rows[K],wp+sz*M,1,r+sz*(inbands_getindex(R.blockstart,K,J) +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
                M=Mpre
                for K=K1+J-BRS1:CS
                    jshft = (j-1)*R.rows[K]  # shift by each column we are to the right
                    BLAS.axpy!(R.rows[K],-2*dt,wp+sz*M,1,r+sz*(inbands_getindex(R.blockstart,K,J) +jshft),1)
                    M += R.rows[K]    # copy all rows in K-th block
                end
            end
        end
    end
    QR.ncols=col
    QR
end
