using ApproxFun

function Jacobi(a,b,N)
  J = zeros(N,N)
  J[1,1]=a[1]
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

function Lmatrix(a,b,N)
  # initial values
  L = zeros(N,N)
  L[1,1] = 1
  L[2,1] = -a[1]/b[1]
  L[2,2] = 1/b[1]
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

n=5
a = rand(n,1)
b = rand(n-1,1)
J = Jacobi(a,b,n)
L = Lmatrix(a,b,2*n)

T=ToeplitzOperator(vec(L[2*n,2*n-1:-1:1]),[L[n,n]])
K = zeros(2*n,2*n)
for i = 1:2*n
  for j = 1:i
    K[i,j] = L[i,j]-T[i,j]
  end
end
K = CompactOperator(K)
