function lyapdiag!(A::Vector,B::Vector,F)
    for i=1:size(F,1),j=1:size(F,2)
        F[i,j] *= 1./(A[i]+B[j])
    end
    
    F
end

# Solve P*Y*T' + S*Y*R' = F 
# where P,T, S and R are upper triangular
function lyapuptriang{N}(P::Array{N},R::Array{N},S::Array{N},T::Array{N},F::Array{N})
    m=size(P,2); n = size(T,2)
    Y=zeros(N,m,n)
    PY = zeros(N,m,m); SY = zeros(N,m,m)
    
    k=n
    while k > 1
        if T[k,k-1] == 0 && R[k,k-1] == 0
            rhs = F[:,k]
            if k < n                
                for j = k+1:n
                    rhs -= R[k,j]*PY[:,j] + T[k,j]*SY[:,j]
                end
            end
    
            Y[:,k] = (R[k,k]*P + T[k,k]*S)\rhs
            PY[:,k] = P*Y[:,k];SY[:,k] = S*Y[:,k]            
            
            k -= 1
        else
            rhs1 = F[:,k-1];rhs2=F[:,k]
            for j = k+1:n
                rhs1 -= R[k-1,j]*PY[:,j] + T[k-1,j]*SY[:,j]
                rhs2 -= R[k,j]*PY[:,j] + T[k,j]*SY[:,j]
            end
            
            SM = zeros(N,2n,2n)
            up=1:n
            down=n+1:2n
            
            SM[up,up]=R[k-1,k-1]*P + T[k-1,k-1]*S
            SM[up,down]=R[k-1,k]*P + T[k-1,k]*S
            SM[down,up] = R[k,k-1]*P + T[k,k-1]*S
            SM[down,down] = R[k,k]*P + T[k,k]*S
            
            UM = SM\[rhs1;rhs2]
            
            Y[:,k-1] = UM[up];Y[:,k]=UM[down]
            
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
    AC=schurfact(A,C)
    BD=schurfact(B,D)
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