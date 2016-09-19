
function Ac_mul_Bpars{RR,T}(A::QROperatorQ{QROperator{RR,Matrix{T},T},T},
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
            resizedata!(A.QR,:,k+M+50)
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


# BLAS apply Q

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