

function Base.nullspace{T}(A::Operator{T};tolerance=10eps(real(T)),maxlength=1_000_000)
   K=transpose_nullspace(qrfact(A'),tolerance,maxlength)
    # drop extra rows, and use QR to determine rank
    Q,R=qr(K,Val{true})
    ind=findfirst(r->abs(r)≤100tolerance,diag(R))
    Kret=ind==0?Q:Q[:,1:ind-1]
    Fun(vec(Kret'),ArraySpace(domainspace(A),1,ind-1))
end


function transpose_nullspace(QR::QROperator,tolerance,maxlength)
    T=eltype(QR)
    resizedata!(QR,:,100)

    m=size(QR.H,1)
    K=zeros(100,m-1)
    K[1:m-1,1:m-1]=eye(m-1)

    for k=1:10
        v=QR.H[:,k]
        QQ=eye(T,m)-2v*v'
        K[:,:]=K*QQ[1:end-1,2:end]
        K[k+m-1,:]=QQ[end,2:end]
    end
    k=11
    α=0.9

    while slnorm(K,floor(Int,k^α),:) >tolerance*k
        if k > maxlength
            warn("Max length of $maxlength reached.")
            break
        end

        if k ≥ QR.ncols
            resizedata!(QR,:,k+100)
        end

        v=QR.H[:,k]
        QQ=(eye(m)-2v*v')
        K[:,:]=K*QQ[1:end-1,2:end]

        if k+m-1 > size(K,1)
            K=pad(K,k+m+100,:)
        end
        K[k+m-1,:]=QQ[end,2:end]
        k+=1
    end

    K[1:floor(Int,k^α),:]
end


Base.nullspace{OO<:Operator}(A::Array{OO};kwds...) = nullspace(interlace(A);kwds...)
