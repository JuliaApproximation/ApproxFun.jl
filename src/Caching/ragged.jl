
CachedOperator(::Type{RaggedMatrix},op::Operator;padding::Bool=false) =
    CachedOperator(op,RaggedMatrix(eltype(op),0,Int[]),padding)

## Grow cached operator

function resizedata!{T<:Number}(B::CachedOperator{T,RaggedMatrix{T}},::Colon,n::Integer)
    if n > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if n > B.datasize[2]
        resize!(B.data.cols,n+1)

        if B.padding
            # K is largest colstop.  We get previous largest by looking at precalulated
            # cols
            K = B.datasize[2]==0?0:B.data.cols[B.datasize[2]+1]-B.data.cols[B.datasize[2]]

            for j = B.datasize[2]+1:n-1
                K = max(K,colstop(B.op,j))
                B.data.cols[j+1] = B.data.cols[j] + K
            end
            K = max(K,colstop(B.op,n))
        else
            for j = B.datasize[2]+1:n-1
                B.data.cols[j+1] = B.data.cols[j] + colstop(B.op,j)
            end
            K = colstop(B.op,n)
        end


        B.data.cols[n+1] = B.data.cols[n] + K
        pad!(B.data.data,B.data.cols[n+1]-1)
        B.data.m = K

        jr=B.datasize[2]+1:n
        kr=1:K
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (K,n)

    end

    B
end

resizedata!{T<:Number}(B::CachedOperator{T,RaggedMatrix{T}},n::Integer,m::Integer) =
    resizedata!(B,:,m)



## Grow QR

QROperator{T}(R::CachedOperator{T,RaggedMatrix{T}}) =
    QROperator(R,RaggedMatrix(T,0,Int[]),0)

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
        for J=1:QR.ncols
            wp=view(W,1:colstop(W,J),J)
            for j=m+1:MO.datasize[2]
                kr=j:j+length(wp)-1
                v=view(MO.data,kr,j)
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


function resizedata!{T<:BlasFloat,MM,DS,RS,BI}(QR::QROperator{CachedOperator{T,RaggedMatrix{T},
                                                                 MM,DS,RS,BI}},
                        ::Colon,col)
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    sz=sizeof(T)

    w=pointer(W.data)
    R=MO.data
    r=pointer(R.data)

    if col ≥ MO.datasize[2]
        m = MO.datasize[2]
        resizedata!(MO,:,col+100)  # double the last rows

        R=MO.data
        r=pointer(R.data)

        # apply previous Householders to new columns of R
        for k=1:QR.ncols
            M=colstop(W,k)  # length of wp
            wp=w+(W.cols[k]-1)*sz  # shift by first index of col J

            for j=m+1:MO.datasize[2]
                v=r+(R.cols[j]+k-2)*sz
                dt=dot(M,wp,1,v,1)
                BLAS.axpy!(M,-2*dt,wp,1,v,1)
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
        w=pointer(W.data)
    end

    for k=QR.ncols+1:col
        cs= colstop(R,k)
        M=cs-k+1

        v=r+sz*(R.cols[k]+k-2)    # diagonal entry of R
        wp=w+sz*(W.cols[k]-1)          # k-th column of W
        BLAS.blascopy!(M,v,1,wp,1)
        W.data[W.cols[k]] += flipsign(BLAS.nrm2(M,wp,1),W.data[W.cols[k]])
        normalize!(M,wp)

        # scale rows entries
        for j=k:MO.datasize[2]
            v=r+(R.cols[j]+k-2)*sz
            dt=dot(M,wp,1,v,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end
    end
    QR.ncols=col
    QR
end





## back substitution

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



## Apply Q


function Ac_mul_Bpars{RR,T}(A::QROperatorQ{QROperator{RR,RaggedMatrix{T},T},T},
                            B::Vector{T},tolerance,maxlength)
    if length(B) > A.QR.ncols
        # upper triangularize extra columns to prepare for \
        resizedata!(A.QR,:,length(B)+size(A.QR.H,1)+10)
    end

    H=A.QR.H
    M=size(H,1)
    m=length(B)
    Y=pad(B,m+M+10)

    k=1
    yp=view(Y,1:length(B))
    while (k ≤ m || norm(yp) > tolerance )
        if k > maxlength
            warn("Maximum length $maxlength reached.")
            break
        end
        if k > A.QR.ncols
            # upper triangularize extra columns to prepare for \
            resizedata!(A.QR,:,k+M+50)
            H=A.QR.H
            M=size(H,1)
        end

        cr=colrange(H,k)

        if k+length(cr)-1>length(Y)
            pad!(Y,2*(k+M))
        end

        wp=view(H,cr,k)
        yp=view(Y,k-1+(cr))

        dt=dot(wp,yp)
        Base.axpy!(-2*dt,wp,yp)
        k+=1
    end
    Fun(resize!(Y,k),domainspace(A))  # chop off zeros
end



# BLAS apply Q

function Ac_mul_Bpars{RR,T<:BlasFloat}(A::QROperatorQ{QROperator{RR,RaggedMatrix{T},T},T},
                            B::Vector{T},tolerance,maxlength)
    if length(B) > A.QR.ncols
        # upper triangularize extra columns to prepare for \
        resizedata!(A.QR,:,length(B)+size(A.QR.H,1)+10)
    end

    H=A.QR.H
    h=pointer(H.data)

    M=size(H,1)
    m=length(B)
    Y=pad(B,m+M+10)

    sz=sizeof(T)

    k=1
    y=pointer(Y)

    yp=y
    while (k ≤ m || BLAS.nrm2(M,yp,1) > tolerance )
        if k > maxlength
            warn("Maximum length $maxlength reached.")
            break
        end
        if k > A.QR.ncols
            # upper triangularize extra columns to prepare for \
            resizedata!(A.QR,:,k+M+50)
            H=A.QR.H
            h=pointer(H.data)
        end


        M=colstop(H,k)

        if k+M-1>length(Y)
            pad!(Y,2*(k+M))
            y=pointer(Y)
        end


        wp=h + sz*(H.cols[k]-1)
        yp=y+sz*(k-1)

        dt=dot(M,wp,1,yp,1)
        BLAS.axpy!(M,-2*dt,wp,1,yp,1)
        k+=1
    end
    Fun(resize!(Y,k),domainspace(A))  # chop off zeros
end
