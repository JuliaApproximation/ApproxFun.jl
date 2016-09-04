#### Caching
function CachedOperator{T}(io::InterlaceOperator{T,1};padding::Bool=false)
    ds=domainspace(io)
    rs=rangespace(io)

    ind=find(op->isinf(size(op,1)),io.ops)
    if length(ind) ≠ 1  || !isbanded(io.ops[ind[1]])  # is almost banded
        return CachedOperator(io,Array(eltype(op),0,0))
    end
    i=ind[1]
    bo=io.ops[i]
    lin,uin=bandwidths(bo)



    # calculate number of rows interlaced
    # each each row in the rangespace increases the lower bandwidth by 1
    nds=0
    md=0
    for k=1:length(io.ops)
        if k ≠ i
            d=dimension(rs[k])
            nds+=d
            md=max(md,d)
        end
    end


    isend=true
    for k=i+1:length(io.ops)
        if dimension(rs[k]) == md
            isend=false
        end
    end

    numoprows=isend?md-1:md
    n=nds+numoprows

    (l,u) = (max(lin+nds,n-1),max(0,uin+1-ind[1]))

    # add extra rows for QR
    if padding
        u+=l
    end


    ret=abzeros(T,n,n+u,l,u,nds)

    # populate the finite rows
    jr=1:n+u
    bcrow=1
    oprow=0
    for k=1:n
        K,J=io.rangeinterlacer[k]

        if K ≠ i
            # fill the fill matrix
            ret.fill.V[:,bcrow] = Matrix(view(io.ops[K],J:J,jr))
            ret.fill.U[k,bcrow] = 1
            bcrow += 1
        else
            oprow+=1
        end


        for j=rowrange(ret.bands,k)
            ret[k,j] = io[k,j]
        end
    end


    CachedOperator(io,ret,(n,n+u),ds,rs,(l,∞))
end


function CachedOperator{T}(io::InterlaceOperator{T,2};padding::Bool=false)
    ds=domainspace(io)
    rs=rangespace(io)
    di=io.domaininterlacer
    ri=io.rangeinterlacer
    ddims=di.iterator.dimensions
    rdims=ri.iterator.dimensions

    # we are only almost banded if every operator is either finite
    # range or banded, and if the # of  ∞ spaces is the same
    # for the domain and range
    isab=all(op->isfinite(size(op,1)) || isbanded(op),io.ops) &&
        count(isinf,ddims) == count(isinf,rdims)

    if !isab
        return CachedOperator(io,Array(T,0,0))
    end


    # these are the bandwidths if we only had ∞-dimensional operators
    d∞=find(isinf,[ddims...])
    r∞=find(isinf,[rdims...])
    p=length(d∞)

    l∞,u∞ = 0,0
    for k=1:p,j=1:p
        l∞=min(l∞,p*bandinds(io.ops[r∞[k],d∞[j]],1)+j-k)
    end
    for k=1:p,j=1:p
        u∞=max(u∞,p*bandinds(io.ops[r∞[k],d∞[j]],2)+j-k)
    end

    # now we move everything by the finite rank
    ncols=mapreduce(d->isfinite(d)?d:0,+,ddims)
    nbcs=mapreduce(d->isfinite(d)?d:0,+,rdims)
    shft=ncols-nbcs
    l∞,u∞=l∞+shft,u∞+shft

    # iterate through finite rows to find worst case bandinds
    l,u=l∞,u∞
        for k=1+nbcs+p,j=1:ncols+p
            N,n=ri[k]
            M,m=di[j]
            l=min(l,bandinds(io.ops[N,M],1)+n-m-k+j)
            u=max(u,bandinds(io.ops[N,M],2)+n-m-k+j)
        end


    # add extra rows for QR
    if padding
        u-=l
    end

    n=1+nbcs+p
    ret=abzeros(T,n,n+u,-l,u,nbcs)

    # populate entries and fill functionals
    bcrow=1
    oprow=0
    jr=1:n+u
    for k=1:n
        K,J=io.rangeinterlacer[k]

        if isfinite(rdims[K] )
            # fill the fill matrix
            ret.fill.V[:,bcrow] = Matrix(view(io,k:k,jr))
            ret.fill.U[k,bcrow] = 1
            bcrow += 1
        else
            oprow+=1
        end


        for j=rowrange(ret.bands,k)
            ret[k,j] = io[k,j]
        end
    end

    CachedOperator(io,ret,(n,n+u),ds,rs,(l,∞))
end

function resizedata!{T<:Number,DS,RS,DI,RI,BI}(co::CachedOperator{T,AlmostBandedMatrix{T},
                                                        InterlaceOperator{T,1,DS,RS,DI,RI,BI}},
                                     n::Integer,::Colon)
    if n ≤ co.datasize[1]
        return co
    end

    (l,u)=bandwidths(co.data.bands)
    pad!(co.data,n,n+u)

    r=rank(co.data.fill)
    ind=findfirst(op->isinf(size(op,1)),co.op.ops)

    k=1
    for (K,J) in co.op.rangeinterlacer
        if K ≠ ind
            co.data.fill.V[co.datasize[2]:end,k] = co.op.ops[K][J,co.datasize[2]:n+u]
            k += 1
            if k > r
                break
            end
        end
    end


    kr=co.datasize[1]+1:n
    jr=max(1,kr[1]-l):n+u
    BLAS.axpy!(1.0,view(co.op.ops[ind],kr-r,jr),
                    view(co.data.bands,kr,jr))

    co.datasize=(n,n+u)
    co
end

function resizedata!{T<:Number,DS,RS,DI,RI,BI}(co::CachedOperator{T,AlmostBandedMatrix{T},
                                                        InterlaceOperator{T,2,DS,RS,DI,RI,BI}},
                                     n::Integer,::Colon)
    if n ≤ co.datasize[1]
        return co
    end

    io=co.op
    ds=domainspace(io)
    rs=rangespace(io)
    di=io.domaininterlacer
    ri=io.rangeinterlacer
    ddims=di.iterator.dimensions
    rdims=ri.iterator.dimensions

    d∞=find(isinf,[ddims...])
    r∞=find(isinf,[rdims...])
    p=length(d∞)


    (l,u)=bandwidths(co.data.bands)
    pad!(co.data,n,n+u)
    co.data
    # r is number of extra rows, ncols is number of extra columns
    r=rank(co.data.fill)
    ncols=mapreduce(d->isfinite(d)?d:0,+,ddims)


    # fill rows
    K=k=1
    while k ≤ r
        if isfinite(dimension(rs[ri[K][1]]))
            co.data.fill.V[co.datasize[2]:end,k] = co.op[K,co.datasize[2]:n+u]
            k += 1
        end
        K += 1
    end

    kr=co.datasize[1]+1:n
    jr=max(1,kr[1]-l):n+u
    io∞=InterlaceOperator(io.ops[r∞,d∞])

    BLAS.axpy!(1.0,view(io∞,kr-r,jr-ncols),view(co.data.bands,kr,jr))

    co.datasize=(n,n+u)
    co
end


resizedata!{T<:Number,DS,RS,DI,RI,BI}(co::CachedOperator{T,AlmostBandedMatrix{T},
                                                        InterlaceOperator{T,1,DS,RS,DI,RI,BI}},
                     n::Integer,m::Integer) = resizedata!(co,max(n,m+bandwidth(co.data.bands,1)),:)


resizedata!{T<:Number,DS,RS,DI,RI,BI}(co::CachedOperator{T,AlmostBandedMatrix{T},
                                                         InterlaceOperator{T,2,DS,RS,DI,RI,BI}},
                      n::Integer,m::Integer) = resizedata!(co,max(n,m+bandwidth(co.data.bands,1)),:)



##
# These give highly optimized routines for delaying with Cached
#

## populate data



function resizedata!{T,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,AlmostBandedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data.bands
    M=R.l+1   # number of diag+subdiagonal bands

    if col+M-1 ≥ MO.datasize[1]
        resizedata!(MO,(col+M-1)+100,:)  # double the last rows
    end

    if col > size(W,2)
        W=QR.H=unsafe_resize!(W,:,2col)
    end

    F=MO.data.fill.U

    for k=QR.ncols+1:col
        W[:,k] = view(R.data,R.u+1:R.u+R.l+1,k) # diagonal and below
        wp=view(W,:,k)
        W[1,k]+= flipsign(norm(wp),W[1,k])
        normalize!(wp)

        # scale banded entries
        for j=k:k+R.u
            dind=R.u+1+k-j
            v=view(R.data,dind:dind+M-1,j)
            dt=dot(wp,v)
            Base.axpy!(-2*dt,wp,v)
        end

        # scale banded/filled entries
        for j=k+R.u+1:k+R.u+M-1
            p=j-k-R.u
            v=view(R.data,1:M-p,j)  # shift down each time
            wp2=view(wp,p+1:M)
            dt=dot(wp2,v)
            for ℓ=k:k+p-1
                @inbounds dt=muladd(conj(W[ℓ-k+1,k]),
                                    unsafe_getindex(MO.data.fill,ℓ,j),dt)
            end
            Base.axpy!(-2*dt,wp2,v)
        end

        # scale filled entries

        for j=1:size(F,2)
            v=view(F,k:k+M-1,j) # the k,jth entry of F
            dt=dot(wp,v)
            Base.axpy!(-2*dt,wp,v)
        end
    end
    QR.ncols=col
    QR
end


function resizedata!{T,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,BandedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data
    M=R.l+1   # number of diag+subdiagonal bands

    if col+M-1 ≥ MO.datasize[1]
        resizedata!(MO,(col+M-1)+100,:)  # double the last rows
    end

    if col > size(W,2)
        W=QR.H=unsafe_resize!(W,:,2col)
    end

    for k=QR.ncols+1:col
        W[:,k] = view(R.data,R.u+1:R.u+R.l+1,k) # diagonal and below
        wp=view(W,:,k)
        W[1,k]+= flipsign(norm(wp),W[1,k])
        normalize!(wp)

        # scale banded entries
        for j=k:k+R.u
            dind=R.u+1+k-j
            v=view(R.data,dind:dind+M-1,j)
            dt=dot(wp,v)
            Base.axpy!(-2*dt,wp,v)
        end

        # scale banded/filled entries
        for j=k+R.u+1:k+R.u+M-1
            p=j-k-R.u
            v=view(R.data,1:M-p,j)  # shift down each time
            wp2=view(wp,p+1:M)
            dt=dot(wp2,v)
            Base.axpy!(-2*dt,wp2,v)
        end
    end
    QR.ncols=col
    QR
end


function resizedata!{T,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,RaggedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    if col ≥ MO.datasize[2]
        m = MO.datasize[2]
        resizedata!(MO,:,col+100)  # double the last rows

        # apply previous Householders to new columns of R
        for J=1:size(W,2)
            wp=view(W,:,J)
            for j=m+1:MO.datasize[2]
                kr=j:j+length(wp)-1
                v=view(MO,kr,j)
                dt=dot(wp,v)
                Base.axpy!(-2*dt,wp,v)
            end
        end
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
        cs = colstop(QR.R,k)
        W[1:cs-k+1,k] = view(MO.data,k:cs,k) # diagonal and below
        wp=view(W,1:cs-k+1,k)
        W[1,k]+= flipsign(norm(wp),W[1,k])
        normalize!(wp)

        # scale rows entries
        kr=k:k+length(wp)-1
        for j=k:MO.datasize[2]
            v=view(MO.data,kr,j)
            dt=dot(wp,v)
            Base.axpy!(-2*dt,wp,v)
        end
    end
    QR.ncols=col
    QR
end


# BLAS versions, requires BlasFloat

function resizedata!{T<:BlasFloat,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,AlmostBandedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data.bands
    M=R.l+1   # number of diag+subdiagonal bands

    if col+M-1 ≥ MO.datasize[1]
        resizedata!(MO,(col+M-1)+100,:)  # double the last rows
    end

    if col > size(W,2)
        W=QR.H=unsafe_resize!(W,:,2col)
    end

    F=MO.data.fill.U

    f=pointer(F)
    m,n=size(R)
    w=pointer(W)
    r=pointer(R.data)
    sz=sizeof(T)
    st=stride(R.data,2)
    stw=stride(W,2)

    for k=QR.ncols+1:col
        v=r+sz*(R.u + (k-1)*st)    # diagonal entry
        wp=w+stw*sz*(k-1)          # k-th column of W
        BLAS.blascopy!(M,v,1,wp,1)
        W[1,k]+= flipsign(BLAS.nrm2(M,wp,1),W[1,k])
        normalize!(M,wp)

        for j=k:k+R.u
            v=r+sz*(R.u + (k-1)*st + (j-k)*(st-1))
            dt=dot(M,wp,1,v,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end

        for j=k+R.u+1:k+R.u+M-1
            p=j-k-R.u
            v=r+sz*((j-1)*st)  # shift down each time
            dt=dot(M-p,wp+p*sz,1,v,1)
            for ℓ=k:k+p-1
                @inbounds dt=muladd(conj(W[ℓ-k+1,k]),
                                    unsafe_getindex(MO.data.fill,ℓ,j),dt)
            end
            BLAS.axpy!(M-p,-2*dt,wp+p*sz,1,v,1)
        end

        fp=f+(k-1)*sz
        fst=stride(F,2)
        for j=1:size(F,2)
            v=fp+fst*(j-1)*sz   # the k,jth entry of F
            dt=dot(M,wp,1,v,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end
    end
    QR.ncols=col
    QR
end



function resizedata!{T<:BlasFloat,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,BandedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    R=MO.data
    M=R.l+1   # number of diag+subdiagonal bands

    if col+M-1 ≥ MO.datasize[1]
        resizedata!(MO,(col+M-1)+100,:)  # double the last rows
    end

    if col > size(W,2)
        W=QR.H=unsafe_resize!(W,:,2col)
    end

    m,n=size(R)
    w=pointer(W)
    r=pointer(R.data)
    sz=sizeof(T)
    st=stride(R.data,2)
    stw=stride(W,2)

    for k=QR.ncols+1:col
        v=r+sz*(R.u + (k-1)*st)    # diagonal entry
        wp=w+stw*sz*(k-1)          # k-th column of W
        BLAS.blascopy!(M,v,1,wp,1)
        W[1,k]+= flipsign(BLAS.nrm2(M,wp,1),W[1,k])
        normalize!(M,wp)

        for j=k:k+R.u
            v=r+sz*(R.u + (k-1)*st + (j-k)*(st-1))
            dt=dot(M,wp,1,v,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end

        for j=k+R.u+1:k+R.u+M-1
            p=j-k-R.u
            v=r+sz*((j-1)*st)  # shift down each time
            dt=dot(M-p,wp+p*sz,1,v,1)
            BLAS.axpy!(M-p,-2*dt,wp+p*sz,1,v,1)
        end
    end
    QR.ncols=col
    QR
end

# back substitution
trtrs!(::Type{Val{'U'}},co::CachedOperator,u::Array) =
                trtrs!(Val{'U'},resizedata!(co,size(u,1),size(u,1)).data,u)

function trtrs!(::Type{Val{'U'}},B::AlmostBandedMatrix,u::Array)
    n=size(u,1)
    A=B.bands
    F=B.fill
    b=bandwidth(A,2)
    nbc = rank(B.fill)
    T=eltype(u)
    pk = zeros(T,nbc)

    for c=1:size(u,2)
        fill!(pk,zero(T))

        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            @simd for j=k+1:n
                @inbounds u[k,c]=muladd(-A.data[k-j+A.u+1,j],u[j,c],u[k,c])
            end

            @inbounds u[k,c] /= A.data[A.u+1,k]
        end

       #filled rows
        for k=n-b-1:-1:1
            @simd for j=1:nbc
                @inbounds pk[j] = muladd(u[k+b+1,c],F.V[k+b+1,j],pk[j])
            end

            @simd for j=k+1:k+b
                @inbounds u[k,c]=muladd(-A.data[k-j+A.u+1,j],u[j,c],u[k,c])
            end

            @simd for j=1:nbc
                @inbounds u[k,c] = muladd(-F.U[k,j],pk[j],u[k,c])
            end

            @inbounds u[k,c] /= A.data[A.u+1,k]
        end
    end
    u
end

function trtrs!(::Type{Val{'U'}},A::BandedMatrix,u::Array)
    n=size(u,1)
    b=bandwidth(A,2)
    T=eltype(u)

    for c=1:size(u,2)
        for k=n:-1:1
            @simd for j=k+1:min(n,k+b)
                @inbounds u[k,c]=muladd(-A.data[k-j+A.u+1,j],u[j,c],u[k,c])
            end

            @inbounds u[k,c] /= A.data[A.u+1,k]
        end
    end
    u
end


function trtrs!(::Type{Val{'U'}},A::RaggedMatrix,u::Array)
    if size(A,1) < size(u,1)
        throw(BoundsError())
    end

    n=size(u,1)
    b=bandwidth(A,2)
    T=eltype(u)

    for c=1:size(u,2)
        for k=n:-1:1
            @inbounds ck = A.cols[k]
            @inbounds u[k,c] /= A.data[ck+k-1]
            BLAS.axpy!(-u[k,c],view(A.data,ck:ck+k-2),view(u,1:k-1,c))
        end
    end
    u
end



## Ac_mul_B! for QROperatorQ

function Ac_mul_Bpars(A::QROperatorQ,B::Vector,tolerance,maxlength)
    T = promote_type(eltype(A),eltype(B))
    Ac_mul_Bpars(Operator{T}(A),Vector{T}(B),tolerance,maxlength)
end


function Ac_mul_Bpars{QR,T}(A::QROperatorQ{QR,T},B::Vector{T},tolerance,maxlength)
    if length(B) > A.QR.ncols
        # upper triangularize extra columns to prepare for \
        resizedata!(A.QR,:,length(B)+size(A.QR.H,1)+10)
    end

    H=A.QR.H
    M=size(H,1)
    m=length(B)
    Y=pad(B,m+M+10)

    k=1
    yp=view(Y,1:M)
    while (k ≤ m+M || norm(yp) > tolerance )
        if k > maxlength
            warn("Maximum length $maxlength reached.")
            break
        end

        if k+M-1>length(Y)
            pad!(Y,2*(k+M))
        end
        if k > A.QR.ncols
            # upper triangularize extra columns to prepare for \
            resizedata!(A.QR,:,2*(k+M))
            H=A.QR.H
        end

        wp=view(H,:,k)
        yp=view(Y,k:k+M-1)

        dt=dot(wp,yp)
        Base.axpy!(-2*dt,wp,yp)
        k+=1
    end
    Fun(resize!(Y,k),domainspace(A))  # chop off zeros
end

# BLAS

function Ac_mul_Bpars{RR,T<:BlasFloat}(A::QROperatorQ{QROperator{RR,Matrix{T},T},T},
                                        B::Vector{T},
                                        tolerance,maxlength)
    if length(B) > A.QR.ncols
        # upper triangularize extra columns to prepare for \
        resizedata!(A.QR,:,length(B)+size(A.QR.H,1)+10)
    end

    H=A.QR.H
    h=pointer(H)

    M=size(H,1)

    b=pointer(B)
    st=stride(H,2)

    sz=sizeof(T)

    m=length(B)
    Y=pad(B,m+M+10)
    y=pointer(Y)

    k=1
    yp=y
    while (k ≤ m+M || BLAS.nrm2(M,yp,1) > tolerance )
        if k > maxlength
            warn("Maximum length $maxlength reached.")
            break
        end

        if k+M-1>length(Y)
            pad!(Y,2*(k+M))
            y=pointer(Y)
        end
        if k > A.QR.ncols
            # upper triangularize extra columns to prepare for \
            resizedata!(A.QR,:,2*(k+M))
            H=A.QR.H
            h=pointer(H)
        end

        wp=h+sz*st*(k-1)
        yp=y+sz*(k-1)

        dt=dot(M,wp,1,yp,1)
        BLAS.axpy!(M,-2*dt,wp,1,yp,1)
        k+=1
    end
    Fun(resize!(Y,k),domainspace(A))  # chop off zeros
end
