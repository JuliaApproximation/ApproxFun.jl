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


backsubstitution!(co::CachedOperator,u::Array) = backsubstitution!(co.data,u)

function backsubstitution!(B::AlmostBandedMatrix,u::Array)
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

function backsubstitution!(A::BandedMatrix,u::Array)
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
