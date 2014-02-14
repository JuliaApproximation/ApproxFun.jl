function lyapdiag!(A::Vector,B::Vector,F)
    for i=1:size(F,1),j=1:size(F,2)
        F[i,j] *= 1./(A[i]+B[j])
    end
    
    F
end

function lyap(A,B,F)
  Λ,V=eig(A)
  Ω,W=eig(B)
    
    F=inv(V)*F*W
    
    Y=lyapdiag!(Λ,Ω,F)
    
    real(V*Y*inv(W))
end