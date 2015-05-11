using ApproxFun

#This now outputs a [-1,1] Chebyshev Fun f such that dmu(s) = f(s)(1-s^2)^{-1/2}
function spectralmeasure(a,b)
  T,K=tkoperators(a,b)
=======
function spectralmeasure(a,b)
    T,K=tkoperators(a,b)
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6

  L = T+K
  Ti=inv(T)
  m=size(K.matrix,1)
  IKin=CompactOperator(inv(full((I+K*Ti)[1:m,1:m]))-eye(m))+I
  Lin=Ti*IKin

<<<<<<< HEAD
  Fun(Lin*[1],Chebyshev)
=======
  x=Fun(identity,[-2,2])
  wΔ=sqrt(4-x^2);wΔ=wΔ/sum(wΔ)
  Fun(Fun(Lin*[1],Ultraspherical{1}([-2,2])),Chebyshev)*wΔ
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6
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


<<<<<<< HEAD
# Finds L such that L(\Delta+Jacobi(a,b-.5))L^{-1} = \Delta, where \Delta = Toeplitz([0,.5])
=======

>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6
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

<<<<<<< HEAD
a = [.2,.3,-.4,.2]
b = [1,1.2,1.3]
T,K = tkoperators(a,b)

=======
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6
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

<<<<<<< HEAD
function jacobioperator(a,b)
    Δ=ToeplitzOperator([.5],[0.,.5])
    B=jacobimatrix(a,b-.5,length(a))
=======
[20. 50.;50. 10.]|>eigvals

α,β=[20.,50.,75.],[1.,1.]
    T,K=tkoperators(α,β)
    λ=eigvals(ApproxFun.companion_matrix(a.coefficients))
    filter!(l->abs(l)≤1.,λ)
    λ=joukowsky(λ)

using AMVW

AMVW.rootsAMVW
Pkg.build("AMVW")
α,β=[20.,50.,75.],[1.,1.]
    T,K=tkoperators(α,β)
    λ=AMVW.rootsAMVW(a.coefficients)
    filter!(l->abs(l)≤1.,λ)
    λ=joukowsky(λ)|>real

sort(λ)
(sort(λ)-λs)./sort(λ)

λ-reverse(λs)

λ

λ=real(joukowsky(complexroots(Fun(T))))

α,β=[20.,50.,75.],[1.,1.]
    T,K=tkoperators(α,β)
    a=Fun(T)
    complexroots(a)


Pkg.add("AMVW")

complexroots(a)
joukowsky(eigvals(ApproxFun.companion_matrix(a.coefficients)))-joukowsky(complexroots(a))

λ[[1,2,3]]-reverse(λs)
T.


reverse(λs)

T,K=tkoperators([0.,0.],[2.])
joukowsky(z)=z+1./z


λ[1]

function jacobioperator(a,b)
    Δ=ToeplitzOperator([1.],[0.,1.])
    B=jacobimatrix(a,b-1,length(a))
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6
    Δ+CompactOperator(B)
end

joukowsky(z)=.5*(z+1./z)

<<<<<<< HEAD
function ql(T::ToeplitzOperator)
  # a(z) = q(z)l(z), where q(z)q(z*) = 1, l is analytic
  # to compute this, take logs: log a(z) = v(z) + ϕ(z), where v(z) + v(z*) = 0, ϕ is analytic
  # then take q = exp(v), l = exp(ϕ)

  # This function is the wrong way around to usual because we want lower triangulars to give analytic symbols
  # In future if T is symmetric then we should be using CosSpace
  a = Fun(ApproxFun.interlace([T.nonnegative[1];T.negative],T.nonnegative[2:end]),Laurent(Circle()))
  la=log(a)
=======
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




>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6

  # The ideal implementation uses SinSpace:
  # v=Fun(2im*la.coefficients[2:2:end],SinSpace(Circle()))

<<<<<<< HEAD
  v=Fun(ApproxFun.interlace([0.;-la.coefficients[2:2:end]],la.coefficients[2:2:end]),Laurent(Circle()))
  q=exp(v)
=======
x=Fun(identity,[-2,2]);w=1/sqrt(4-x^2)
    lanczos(w,10)
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6

  Qtranspose = ToeplitzOperator(q)
  ToeplitzOperator(Qtranspose.nonnegative[2:end],[Qtranspose.nonnegative[1];Qtranspose.negative])

<<<<<<< HEAD
  ϕ=la-v
  l=exp(ϕ)
  Ltranspose = ToeplitzOperator(l)
  ToeplitzOperator(Ltranspose.nonnegative[2:end],[Ltranspose.nonnegative[1]])

end

T = ToeplitzOperator([2.,.5])
T[1:10,1:10]
a = Fun(ApproxFun.interlace([T.nonnegative[1];T.negative],T.nonnegative[2:end]),Laurent(Circle()))
  la=log(a)
v=Fun(ApproxFun.interlace([0.;-la.coefficients[2:2:end]],la.coefficients[2:2:end]),Laurent(Circle()))
  q=exp(v)
Qtranspose = ToeplitzOperator(q)
  Q = ToeplitzOperator(Qtranspose.nonnegative[2:end],[Qtranspose.nonnegative[1];Qtranspose.negative])

(Qtranspose*Q-I)[2:20,2:20]|>full|>norm
(Q*Qtranspose-I)[9:100,9:100]|>full|>norm
=======
wL=Fun(1,[-2,2]);wL=wL/sum(wL)

@code_typed spectralmeasure(zeros(n),lanczos(wL,n-1)[2])

m=length(w)

n=60; @time w=spectralmeasure(0.01(rand(n)-0.5),ones(n-1))
    ApproxFun.plot(w)

(lanczos(w,n-1)[2]-lanczos(wL,n-1)[2])|>norm


length(w)
lanczos(w,3)
lanczos(w,3)[2]-[b;1]






######
# QL For Toeplitz
######

using ApproxFun

a=Fun(ζ->ζ+1/ζ+2.7,Laurent(Circle()))

ToeplitzOperator(a)[1:10,1:10]

la=log(a)

ζ=exp(.1im);la[ζ]+la[1/ζ]

v=Fun(2im*la.coefficients[2:2:end],SinSpace(Circle()))

v=Fun(ApproxFun.interlace([0.;-la.coefficients[2:2:end]],la.coefficients[2:2:end]),
        Laurent(Circle()))
v[ζ]+v[1/ζ]

φ=la-v
l=exp(φ)
q=exp(v)

q*Φ-a|>norm

q[ζ]*q[1/ζ]




A=ToeplitzOperator(a)'
Q=ToeplitzOperator(q)'
L=ToeplitzOperator(l)'

using SO
L[1:10,1:10]|>chopm

q[ζ]*q[1/ζ]

(Q*(Q'))[1:30,1:30]|>chopm

QM[:,5]|>chopm

QM=Q[1:50,1:50]
QM*QM'|>chopm

(A-Q*L)[1:10,1:10]|>chopm

(ToeplitzOperator(q)*ToeplitzOperator(l)

Φ.coefficients
>>>>>>> 8ab3153624e84aa12fbc92ff35524f9663ed61f6
