
normalize!(w) = BLAS.scal!(length(w),inv(norm(w)),w,1)
function normalize!(n,w)
    BLAS.scal!(n,inv(BLAS.nrm2(n,w,1)),w,1)
end

function resizecols!(W::Matrix,m)
    n=size(W,1)
    reshape(resize!(vec(W),n*m),n,m)
end


type QROperator{B,S,T}  <: Operator{T}
    R::MutableOperator{T,B,S}
    H::Matrix{T} # Contains the Householder reflections
    ncols::Int   # number of cols already upper triangularized
end


immutable QROperatorR{QRT,T} <: Operator{T}
    QR::QRT
end

QROperatorR(QR) = QROperatorR{typeof(QR),eltype(QR)}(QR)


immutable QROperatorQ{QRT,T} <: Operator{T}
    QR::QRT
end

QROperatorQ(QR) = QROperatorQ{typeof(QR),eltype(QR)}(QR)

function getindex(QR::QROperator,d::Symbol)
    d==:Q && return QROperatorQ(QR)
    d==:R && return QROperatorR(QR)

    error("Symbol not recognized")
end


function Base.qrfact{OO<:Operator}(A::Vector{OO})
    R = MutableOperator(A)
    M=R.data.l+1   # number of diag+subdiagonal bands
    H=Array(mapreduce(eltype,promote_type,A),M,100)
    QROperator(R,H,0)
end

function Base.qr{OO<:Operator}(A::Vector{OO})
    QR=qrfact(A)
    QR[:Q],QR[:R]
end



## populate data


function resizedata!(QR::QROperator,col)
    MO=QR.R
    W=QR.H

    T=eltype(QR)
    M=MO.data.l+1   # number of diag+subdiagonal bands

    if col+M-1 ≥ MO.datalength
        resizedata!(MO,(col+M-1)+100)  # double the last rows
    end

    if col > size(W,2)
        W=QR.H=resizecols!(W,2col)
    end


    R=MO.data
    F=MO.fill.data

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
        W[1,k]+= sign(W[1,k])*BLAS.nrm2(M,wp,1)
        normalize!(M,wp)

        for j=k:k+R.u
            v=r+sz*(R.u + (k-1)*st + (j-k)*(st-1))
            dt=BLAS.dot(M,v,1,wp,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end

        for j=k+R.u+1:k+R.u+M-1
            p=j-k-R.u
            v=r+sz*((j-1)*st)  # shift down each time
            dt=BLAS.dot(M-p,v,1,wp+p*sz,1)
            for ℓ=k:k+p-1
                @inbounds dt+=conj(W[ℓ-k+1,k])*unsafe_getindex(MO.fill,ℓ,j)
            end
            BLAS.axpy!(M-p,-2*dt,wp+p*sz,1,v,1)
        end

        fp=f+(k-1)*sz
        fst=stride(F,2)
        for j=1:size(F,2)
            v=fp+fst*(j-1)*sz   # the k,jth entry of F
            dt=BLAS.dot(M,v,1,wp,1)
            BLAS.axpy!(M,-2*dt,wp,1,v,1)
        end
    end
    QR.ncols=col
    QR
end




## Multiplication routines


Base.At_mul_B{T<:Real}(A::QROperatorQ{T},B::Union{Vector{T},Matrix{T}}) = Ac_mul_B(A,B)
Base.At_mul_B!{T<:Real}(Y,A::QROperatorQ{T},B::Union{Vector{T},Matrix{T}}) = Ac_mul_B!(Y,A,B)


# Each Householder increases the size og B

function Base.Ac_mul_B{QR,T<:BlasFloat}(A::QROperatorQ{QR,T},B::Vector{T})
    if length(B) > A.QR.ncols
        # upper triangularize extra columns to prepare for \
        resizedata!(A.QR,length(B)+size(A.QR.H,1)+10)
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
    tol=eps()/10

    k=1
    yp=y
    while k ≤ m+M || BLAS.nrm2(M,yp,1) > tol
        if k+M-1>length(Y)
            pad!(Y,2*(k+M))
            y=pointer(Y)
        end
        if k > A.QR.ncols
            # upper triangularize extra columns to prepare for \
            resizedata!(A.QR,2*(k+M))
            H=A.QR.H
            h=pointer(H)
        end

        wp=h+sz*st*(k-1)
        yp=y+sz*(k-1)

        dt=BLAS.dot(M,yp,1,wp,1)
        BLAS.axpy!(M,-2*dt,wp,1,yp,1)
        k+=1
    end
    Y
end


# function Base.\(R::QROperatorR,b::Vector)
#     if length(b) > R.QR.ncols
#         # upper triangularize columns
#         resizedata!(R,length(b))
#     end
# s
# end
