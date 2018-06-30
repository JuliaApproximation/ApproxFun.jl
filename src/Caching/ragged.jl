
CachedOperator(::Type{RaggedMatrix},op::Operator;padding::Bool=false) =
    CachedOperator(op,RaggedMatrix{eltype(op)}(undef, 0, Int[]),padding)

## Grow cached operator

function resizedata!(B::CachedOperator{T,RaggedMatrix{T}},::Colon,n::Integer) where T<:Number
    if n > size(B,2)
        throw(ArgumentError("Cannot resize beyond size of operator"))
    end

    if n > B.datasize[2]
        resize!(B.data.cols,n+1)

        if B.padding
            # K is largest colstop.  We get previous largest by looking at precalulated
            # cols
            K = B.datasize[2]==0 ? 0 : B.data.cols[B.datasize[2]+1]-B.data.cols[B.datasize[2]]

            for j = B.datasize[2]+1:n
                K = max(K,colstop(B.op,j))
                B.data.cols[j+1] = B.data.cols[j] + K
            end
        else
            K = B.datasize[2]==0 ? 0 : B.data.m# more robust but slower: maximum(diff(B.data.cols))

            for j = B.datasize[2]+1:n
                cs = colstop(B.op,j)
                K = max(K,cs)
                B.data.cols[j+1] = B.data.cols[j] + cs
            end
        end

        # avoid padding with negative length
        if B.data.cols[n+1] ≤ 0
            return B
        end

        pad!(B.data.data,B.data.cols[n+1]-1)
        B.data.m = K

        jr=B.datasize[2]+1:n
        kr=1:K
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (K,n)

    end

    B
end

function resizedata!(B::CachedOperator{T,RaggedMatrix{T}},n::Integer,m::Integer) where T<:Number
    resizedata!(B,:,m)
    B.data.m = max(B.data.m,n)   # make sure we have at least n rows

    B
end


## Grow QR

QROperator(R::CachedOperator{T,RaggedMatrix{T}}) where {T} =
    QROperator(R,RaggedMatrix{T}(undef,0,Int[]),0)

function resizedata!(QR::QROperator{CachedOperator{T,RaggedMatrix{T},
                                                  MM,DS,RS,BI}},
         ::Colon,col) where {T,MM,DS,RS,BI}
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    if col > MO.datasize[2]
        m = MO.datasize[2]
        resizedata!(MO,:,col+100)  # last rows plus a bunch more

        # apply previous Householders to new columns of R
        for J=1:QR.ncols
            wp=view(W,1:colstop(W,J),J)
            for j=m+1:MO.datasize[2]
                kr=J:J+length(wp)-1
                v=view(MO.data,kr,j)
                dt=dot(wp,v)
                LinearAlgebra.axpy!(-2*dt,wp,v)
            end
        end
    end


    if col > size(W,2)
        m=size(W,2)
        resize!(W.cols,col+101)

        for j=m+1:col+100
            cs=colstop(MO.data,j)
            W.cols[j+1]=W.cols[j] + cs-j+1
            W.m=max(W.m,cs-j+1)
        end

        resize!(W.data,W.cols[end]-1)
    end

    for k=QR.ncols+1:col
        cs = colstop(MO.data,k)
        W[1:cs-k+1,k] = view(MO.data,k:cs,k) # diagonal and below
        wp=view(W,1:cs-k+1,k)
        W[1,k]+= flipsign(norm(wp),W[1,k])
        normalize!(wp)

        # scale rows entries
        kr=k:k+length(wp)-1
        for j=k:MO.datasize[2]
            v=view(MO.data,kr,j)
            dt=dot(wp,v)
            LinearAlgebra.axpy!(-2*dt,wp,v)
        end
    end
    QR.ncols=col
    QR
end


# BLAS versions, requires BlasFloat


function resizedata!(QR::QROperator{CachedOperator{T,RaggedMatrix{T},
                                       MM,DS,RS,BI}},
::Colon,col) where {T<:BlasFloat,MM,DS,RS,BI}
    if col ≤ QR.ncols
        return QR
    end

    MO=QR.R
    W=QR.H

    sz=sizeof(T)

    w=pointer(W.data)
    R=MO.data
    r=pointer(R.data)

    if col > MO.datasize[2]
        m = MO.datasize[2]
        resizedata!(MO,:,col+100)  # last rows plus a bunch more

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
        resize!(W.cols,col+101)

        for j=m+1:col+100
            cs=colstop(R,j)
            q_len=cs-j+1  # number of entries in j-th column manipulated
            @assert q_len > 0   # Otherwise, diagonal is not included
            W.cols[j+1]=W.cols[j] + q_len
            W.m=max(W.m,q_len)
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
for ArrTyp in (:AbstractVector, :AbstractMatrix)
    @eval function ldiv!(U::UpperTriangular{T, SubArray{T, 2, RaggedMatrix{T}, Tuple{UnitRange{Int}, UnitRange{Int}}, false}},
                             u::$ArrTyp{T}) where T
        n = size(u,1)
        n == size(U,1) || throw(DimensionMismatch())

        V = parent(U)
        @assert parentindexes(V)[1][1] == 1
        @assert parentindexes(V)[2][1] == 1

        A = parent(V)

        for c=1:size(u,2)
            for k=n:-1:1
                @inbounds ck = A.cols[k]
                @inbounds u[k,c] /= A.data[ck+k-1]
                BLAS.axpy!(-u[k,c], view(A.data,ck:ck+k-2), view(u,1:k-1,c))
            end
        end
        u
    end
end



## Apply Q


function mulpars(Ac::Adjoint{T,<:QROperatorQ{QROperator{RR,RaggedMatrix{T},T},T}},
                      B::AbstractVector{T},tolerance,maxlength) where {RR,T}
    A = parent(Ac)
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
    while (k ≤ m || norm(yp) > tolerance )
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
        LinearAlgebra.axpy!(-2*dt,wp,yp)
        k+=1
    end
    resize!(Y,k)  # chop off zeros
end



# BLAS apply Q

function mulpars(Ac::Adjoint{T,<:QROperatorQ{QROperator{RR,RaggedMatrix{T},T},T}},
           B::AbstractVector{T},tolerance,maxlength) where {RR,T<:BlasFloat}
    A = parent(Ac)
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
    while (k ≤ m || BLAS.nrm2(M,yp,1) > tolerance )
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
    resize!(Y,k)  # chop off zeros
end
