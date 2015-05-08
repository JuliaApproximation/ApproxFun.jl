using ApproxFun

function spectralmeasure(a,b)
    T,K=tkoperators(a,b)

  L = T+K
  Ti=inv(T)
  m=size(K.matrix,1)
  IKin=CompactOperator(inv(full((I+K*Ti)[1:m,1:m]))-eye(m))+I
  Lin=Ti*IKin

  x=Fun(identity,[-2,2])
  wΔ=sqrt(4-x^2);wΔ=wΔ/sum(wΔ)
  Fun(Fun(Lin*[1],Ultraspherical{1}([-2,2])),Chebyshev)*wΔ
end

function tkoperators(a,b)
    assert(length(a)-length(b)==1)
    n = length(a)
    L = Lmatrix(a,b,2n)

    T=ToeplitzOperator(vec(L[2*n,2*n-1:-1:1]),[L[n,n]])
    K = zeros(2*n,2*n)
    for i = 1:2*n
        for j = 1:i
            K[i,j] = L[i,j]-T[i,j]
        end
    end

    K = CompactOperator(K)
    T,K
end



function Lmatrix(a,b,N)
  # initial values
  L = zeros(N,N)
  L[1,1] = 1
  L[2,1] = -a[1]/b[1]
  L[2,2] = 1/b[1]

  n = length(a)
  # the generic case.
  for i = 3:n
      L[i,1] = (L[i-1,2]-a[i-1]*L[i-1,1]-b[i-2]*L[i-2,1])/b[i-1]
      for j = 2:i
          L[i,j] = (L[i-1,j+1]-a[i-1]*L[i-1,j]+L[i-1,j-1]-b[i-2]*L[i-2,j])/b[i-1]
      end
  end
  # this special case happens because b[n] = 1 but a[n] != 0
  L[n+1,1] = L[n,2]-a[n]*L[n,1]-b[n-1]*L[n-1,1]
  for j = 2:n+1
      L[n+1,j] = L[n,j+1]-a[n]*L[n,j]+L[n,j-1]-b[n-1]*L[n-1,j]
  end
  # this case is where b[m],b[m-1],b[m-2] = 1, and a[m],a[m-1] = 0, like Chebyshev
  for m = n+2:N
      L[m,1] = L[m-1,2]-L[m-2,1]
      for j = 2:m-1
        L[m,j] = L[m-1,j+1]-L[m-2,j]+L[m-1,j-1]
      end
      L[m,m] = -L[m-2,m] + L[m-1,m-1]
  end
  L
end

function jacobimatrix(a,b,N)
    J = zeros(N,N)
    J[1,1]=a[1]
    n=length(a)
    for i = 1:n-1
        J[i+1,i+1] = a[i+1]
        J[i,i+1] = b[i]
        J[i+1,i] = b[i]
    end
    for i = n:N-1
        J[i,i+1] = 1
        J[i+1,i] = 1
    end
    J
end


T,K=tkoperators([0.,0.],[2.])
joukowsky(z)=z+1./z
λ=real(joukowsky(complexroots(Fun(T))))

λ[1]

function jacobioperator(a,b)
    Δ=ToeplitzOperator([1.],[0.,1.])
    B=jacobimatrix(a,b-1,length(a))
    Δ+CompactOperator(B)
end


a,b=[2.,0.],[1.]
    J=jacobioperator(a,b)
    T,K=tkoperators(a,b)


#####
# QL Decomposition of J around λ[1] = 2.5

λ=joukowsky(complexroots(Fun(T)))|>real
n=20;    A=J[1:n,1:n]|>full
    L=A-λ[1]*I|>qlgivens!
    Q=(A-λ[1]*I)*inv(L)
    Q*Q'
L*A*inv(L)|>chopm

L=ToeplitzOperator([-2.,0.5],[2.])
T=ToeplitzOperator([1],[-λ[1],1.])


(T*inv(L))
 q
cfs=(T*inv(L))[1:148,2]


q=Fun([cfs[2:end]],Taylor())⊕Fun([cfs[1]],Hardy{false}())
qb=Fun(z->q[1/z],space(q))
space(q)

q*qb


ToeplitzOperator(q)[1:10,1:10]

cfs[1]

q.coefficients

Fun(T)

Fun(inv(L)).coefficients

(J-λ[1]*I)[1:10,1:10]

L

w=spectralmeasure([0.5,0.],[1.])

ApproxFun.plot(w)

J=jacobioperator([0.,0.],[2.])
n=200;    A=J[1:n,1:n]|>full
    L=A|>copy|>qlgivens!
    L*A*inv(L)|>chopm

Δ=ToeplitzOperator([1.],[0.,1.])
    B=jacobimatrix([0.,0.],[2.-1.],2)
    J=Δ+CompactOperator(B)
    J[1:10,1:10]


jacobimatrix([0.,0.],[2.],100)|>eigvals
T[1:10,1:10]|>full|>inv






x=Fun(identity,[-2,2]);w=1/sqrt(4-x^2)
    lanczos(w,10)


wL=Fun(1,[-2,2]);wL=wL/sum(wL)

@code_typed spectralmeasure(zeros(n),lanczos(wL,n-1)[2])

m=length(w)

n=60; @time w=spectralmeasure(0.01(rand(n)-0.5),ones(n-1))
    ApproxFun.plot(w)

(lanczos(w,n-1)[2]-lanczos(wL,n-1)[2])|>norm


length(w)
lanczos(w,3)
lanczos(w,3)[2]-[b;1]


