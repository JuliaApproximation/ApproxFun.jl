using ApproxFun

#This now outputs a [-1,1] Chebyshev Fun f that is weighted by (1-x^2)^{-.5} such that dmu(s) = f(s)
function spectralmeasure(a,b)
  T,K=tkoperators(a,b)
  L = T+K
  Ti=inv(T)
  m=size(K.matrix,1)
  IKin=CompactOperator(inv(full((I+K*Ti)[1:m,1:m]))-eye(m))+I
  Lin=Ti*IKin
  Fun(Lin,JacobiWeight(-0.5,-0.5,Chebyshev()))
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

# Finds L such that L(\Delta+Jacobi(a,b-.5))L^{-1} = \Delta, where \Delta = Toeplitz([0,.5])
function Lmatrix(a,b,N)
  # initial values
  L = zeros(N,N)
  L[1,1] = 1
  L[2,1] = -a[1]/b[1]
  L[2,2] = 1/b[1]

  n = length(a)
  # the generic case.
  for i = 3:n
      L[i,1] = (L[i-1,2]/2-a[i-1]*L[i-1,1]-b[i-2]*L[i-2,1])/b[i-1]
      for j = 2:i
          L[i,j] = (L[i-1,j+1]/2-a[i-1]*L[i-1,j]+L[i-1,j-1]/2-b[i-2]*L[i-2,j])/b[i-1]
      end
  end
  # this special case happens because b[n] = 1/2 but a[n] != 0
  L[n+1,1] = L[n,2]-2*a[n]*L[n,1]-2*b[n-1]*L[n-1,1]
  for j = 2:n+1
      L[n+1,j] = L[n,j+1]-2*a[n]*L[n,j]+L[n,j-1]-2*b[n-1]*L[n-1,j]
  end
  # this case is where b[m],b[m-1],b[m-2] = 1/2, and a[m],a[m-1] = 0, like Chebyshev
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
        J[i,i+1] = .5
        J[i+1,i] = .5
    end
    J
end


function jacobioperator(a,b)
    Δ=ToeplitzOperator([.5],[0.,.5])
    B=jacobimatrix(a,b-.5,length(a))
    Δ+CompactOperator(B)
end

joukowsky(z)=.5*(z+1./z)


function givenstail(t0,t1)
    s∞ = (t0 - sqrt(t0^2-4t1^2))/(2t1)
    c∞ = -sqrt(1-s∞^2)
    α = t1*c∞
    β = c∞*t0 - s∞*α
    l0 = (t0 + sqrt(t0^2-4t1^2))/(2t1)
    l1 = 2t1
    l2 = t1*s∞
    s∞,c∞,α,β,ToeplitzOperator([l1,l2],[l0])
end


function tridql!(L::Matrix)
    n=size(L,1)

  # Now we do QL for the compact part in the top left
    Q = eye(eltype(L),n)
    for i = n:-1:2
        nrm=sqrt(L[i-1,i]^2+L[i,i]^2)
        c,s = L[i,i]/nrm, L[i-1,i]/nrm
        if i > 2
            L[i-1:i,i-2:i] = [c -s; s c]*L[i-1:i,i-2:i]
            L[i-1,i]=0
        else
            L[i-1:i,i-1:i] = [c -s; s c]*L[i-1:i,i-1:i]
            L[i-1,i]=0
        end
        G=eye(eltype(L),n)
        G[i-1:i,i-1:i]=[c s; -s c]
        Q = Q*G
    end
    Q,L
end

import ApproxFun:BandedOperator,ToeplitzOperator

immutable ToeplitzGivens <: BandedOperator{Float64}
    c::Float64
    s::Float64
end
s
bandinds(T::ToeplitzGivens)=-ceil(Int,(-36-2log(abs(c)))/log(abs(s))),1

function ToeplitzOperator(T::ToeplitzGivens)
    c,s=T.c,T.s
    nonneg=[c^2,s]
    m=-bandinds(T,1)
    if m ≥ 1
        neg=Array(Float64,m)
        neg[1]=-s*nonneg[1]
        for k=2:m
           neg[k]=-s*neg[k-1]
        end
    else
        neg=[]
    end
    ToeplitzOperator(neg,nonneg)
end


t0=2.;t1=0.5;
s,c,α,β,TL=givenstail(t0,t1)

c,s

c^2
s,c^2,-s*c^2,s^2*c^2,(-s)^3*c^2

*s^k ≤ ε/c^2
k*log(s) ≤ log(ε)-2log(c)
k ≥  (log(ε)-2log(abs(c)))/log(abs(s))

k=ceil(Int,(-36-2log(abs(c)))/log(abs(s)))

c

s^k

log(eps())



function ql(a,b,t0,t1,N)
    @assert t0^2>=4t1^2
    # The Givens rotations coming from infinity (with parameters c∞ and s∞) leave us with the almost triangular
    # a[n-1]  b[n-1]   0    0    0
    # b[n-1]   a[n]   t1    0    0
    #   0       α      β    0    0
    #   0      l2     l1   l0    0
    #   0       0     l2   l1   l0

    s∞,c∞,α,β,TL=givenstail(t0,t1)

    # Here we construct this matrix as L
    n = length(a)
    L = jacobimatrix(a,b,n+1)
    L[n,n+1] = t1
#    L[n+1,n+2] = 0
    L[n+1,n+1]=β
    L[n+1,n]=α

    Q,L=tridql!(L)

    Q2=eye(N)
    Q2[1:n+1,1:n+1]=Q
    Q=Q2


    # Now we add the Givens rotations at infinity to Qt
    for i = n+2:N
        G=eye(N)
        G[i-1:i,i-1:i]=[c∞ s∞; -s∞ c∞]
        Q = G*Q
    end
    Q,L
end

ql([1.,2.,3.],[5.,2.,1.],2.,0.5,10)




# function ql(T::ToeplitzOperator)
#  # a(z) = q(z)l(z), where q(z)q(z*) = 1, l is analytic
#  # to compute this, take logs: log a(z) = v(z) + ϕ(z), where v(z) + v(z*) = 0, ϕ is analytic
#  # then take q = exp(v), l = exp(ϕ)
#
#  a = Fun(ApproxFun.interlace([T.nonnegative[1];T.negative],T.nonnegative[2:end]),Laurent(Circle()))
#  la=log(a)
#
#  v=Fun(ApproxFun.interlace([0.;-la.coefficients[2:2:end]],la.coefficients[2:2:end]),Laurent(Circle()))
#  q=exp(v)
#  Qtranspose = ToeplitzOperator(q)
#  Q = ToeplitzOperator(Qtranspose.nonnegative[2:end],[Qtranspose.nonnegative[1];Qtranspose.negative])
#
#  ϕ=la-v
#  l=exp(ϕ)
#  Ltranspose = ToeplitzOperator(l)
#  L = ToeplitzOperator(Ltranspose.nonnegative[2:end],[Ltranspose.nonnegative[1]])
#  Q,L
# end
