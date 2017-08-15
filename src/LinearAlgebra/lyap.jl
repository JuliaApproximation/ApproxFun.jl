function lyapdiag!(A::Vector,B::Vector,F)
    for i=1:size(F,1),j=1:size(F,2)
        F[i,j] *= 1./(A[i]+B[j])
    end

    F
end

# Solve P*Y*R' + S*Y*T' = F
# where P,R,S and T are upper triangular
function lyapuptriang(P,R,S,T,F::AbstractArray{N}) where N
    m=size(P,2); n = size(T,2)
    Y=Matrix{N}(m,n)
    PY = Matrix{N}(m,n); SY = Matrix{N}(m,n)

    k=n
    while k > 1
        if T[k,k-1] == 0 && R[k,k-1] == 0
            rhs = F[:,k]
            if k < n
                for j = k+1:n, l=1:m
                    rhs[l] -= R[k,j]*PY[l,j] + T[k,j]*SY[l,j]
                end
            end

            Y[:,k] = (R[k,k]*P + T[k,k]*S)\rhs
            PY[:,k] = P*Y[:,k];SY[:,k] = S*Y[:,k]

            k -= 1
        else
            rhs1 = F[:,k-1];rhs2=F[:,k]
            for j = k+1:n, l=1:m
                rhs1[l] -= R[k-1,j]*PY[l,j] + T[k-1,j]*SY[l,j]
                rhs2[l] -= R[k,j]*PY[l,j] + T[k,j]*SY[l,j]
            end

            SM = Matrix{N}(2m,2m)

            for i=1:m,j=1:m
                SM[i,j] = R[k-1,k-1]*P[i,j] + T[k-1,k-1]*S[i,j]
                SM[i,j+m] = R[k-1,k]*P[i,j] + T[k-1,k]*S[i,j]
                SM[i+m,j] = R[k,k-1]*P[i,j] + T[k,k-1]*S[i,j]
                SM[i+m,j+m] = R[k,k]*P[i,j] + T[k,k]*S[i,j]
            end

            UM = SM\[rhs1;rhs2]

            Y[:,k-1] = UM[1:m];Y[:,k]=UM[m+1:end]

            PY[:,k]=P*Y[:,k];PY[:,k-1]=P*Y[:,k-1]
            SY[:,k]=S*Y[:,k];SY[:,k-1]=S*Y[:,k-1]

            k -=2
        end
    end

    if k == 1
        rhs =F[:,1]
        PY[:,2]=P*Y[:,2]
        SY[:,2]=S*Y[:,2]
        for j=2:n
            rhs -= R[1,j]*PY[:,j] + T[1,j]*SY[:,j]
        end

        Y[:,1]=(R[1,1]*P + T[1,1]*S)\rhs
    end
    Y
end


##Solves A*X*B.' + C*X*D.' = E
function lyap(A,B,C,D,E)
    AC=schurfact(full(A),full(C))
    BD=schurfact(full(B),full(D))
    Q1=AC[:left];Q2=BD[:left]
    Z1=AC[:right];Z2=BD[:right]

    Y=lyapuptriang(AC[:S],BD[:S],AC[:T],BD[:T],Q1'*E*Q2)
    Z1*Y*Z2'
end


# function lyap(A,B,F)
#   Λ,V=eig(A)
#   Ω,W=eig(B)
#
#     F=inv(V)*F*W
#
#     Y=lyapdiag!(Λ,Ω,F)
#
#     real(V*Y*inv(W))
# end
