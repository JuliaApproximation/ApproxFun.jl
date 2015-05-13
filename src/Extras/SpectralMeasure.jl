using ApproxFun
    import ApproxFun:BandedOperator,ToeplitzOperator,tridql!,bandinds,DiracSpace

joukowsky(z)=.5*(z+1./z)

function spectralmeasure(a,b)
    T,K=tkoperators(a,b)
    Tfun = Fun([T[1,1];T.negative],Taylor)
    eigs=real(map!(joukowsky,filter!(z->abs(z)<1,complexroots(Tfun))))
    @assert length(eigs) ≤ 1 #TODO: fix
    if isempty(eigs)
        L=T+K
        coeffs = L\[1]
        Fun(coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))
    else
        t1,t0=0.5,-maximum(eigs)

        Q,L1=ql(a+t0,b,t0,t1)
        LQ=L1*Q
        acont=Float64[LQ[k,k] for k=1:length(a)+1]-t0;
        bcont=Float64[LQ[k,k+1] for k=1:length(a)];
        L2=contspectraltransform(acont[2:end],bcont[2:end])

        q0=Q'*[1]
        coeffs=L2\q0[2:end]
        μ1=Fun(coeffs,JacobiWeight(0.5,0.5,Ultraspherical{1}()))
        c=q0[2]
        Fun([q0[1]^2;coeffs],DiracSpace(JacobiWeight(0.5,0.5,Ultraspherical{1}()),[-t0]))
    end
end

function contspectraltransform(a,b)
    T,K=tkoperators(a,b)
    T+K
end

function contspectralmeasure(a,b)
    L=contspectraltransform(a,b)
    coeffs = L\[1]
    Fun(coeffs,JacobiWeight(.5,.5,Ultraspherical{1}())),L
end



α=1.91;a=[-.032145643-α,0.1];b=[0.5];
    spectralmeasure(a,b)|>ApproxFun.plot



μ,L2=contspectralmeasure(acont[2:end],bcont[2:end])

J0=jacobioperator(a,b,0.,0.5)



L2'*q0[2:end]

q0[1]^2+(c*Fun(μ1,JacobiWeight(0.5,0.5,Chebyshev()))*2/π|>sum)


μ=Fun([q0[1]^2;coeffs],DiracSpace(JacobiWeight(0.5,0.5,Ultraspherical{1}()),[-t0]))

f=Fun(exp,Ultraspherical{2}())
sum(f)
sum(μ1)
ApproxFun.plot(μ)
ApproxFun.plot(μ)




q0[2:end]


μ2

ApproxFun.plot(μ2)

inv(L1)

t0,t1

Q[1:10,1:10]
(Q*L1-jacobioperator(a+t0,b,t0,0.5))[1:10,1:10

jacobioperator(a,b,0.,0.5)[1:100,1:100]|>full|>eigvals

(LQ[1:100,1:100]|>full|>eigvals)-t0


L1[1:10,1:10]

(Q*L1)[1:10,1:10]

L1[1:10,1:10]

L1*[1]

ApproxFun.plot(μ)

acont

(L*Q)[1:10,1:10]|>chopm

function tkoperators(a,b)
    assert(length(a)-length(b)==1)
    n = length(a)
    L = Lmatrix(a,b,2n)

    T=ToeplitzOperator(vec(L[2*n,2*n-1:-1:1]),[L[2*n,2*n]])
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
    n = length(a)
    @assert n-length(b)==1
    bext = [b; .5]
    L = zeros(N,N)
    L[1,1] = 1
    L[2,1] = -a[1]/bext[1]
    L[2,2] = 0.5/bext[1]
    # for Chebyshev T, use L[2,2] = 1/(sqrt(2)*bext[1])

    # the generic case.
    for i = 3:n+1
        L[i,1] = (L[i-1,2]/2-a[i-1]*L[i-1,1]-bext[i-2]*L[i-2,1])/bext[i-1]
        # for Chebyshev T, use
        # L[i,1] = (L[i-1,2]/sqrt(2)-a[i-1]*L[i-1,1]-bext[i-2]*L[i-2,1])/bext[i-1]
        # L[i,2] = (L[i-1,3]/2+L[i-1,1]/sqrt(2)-a[i-1]*L[i-1,2]-bext[i-2]*L[i-2,2])/bext[i-1]
        # and go from j=3:i
        for j = 2:i
            L[i,j] = (L[i-1,j+1]/2+L[i-1,j-1]/2-a[i-1]*L[i-1,j]-bext[i-2]*L[i-2,j])/bext[i-1]
        end
    end
    # this case is where b[m],b[m-1],b[m-2] = 1/2, and a[m],a[m-1] = 0, like Chebyshev
    for m = n+2:N
        L[m,1] = L[m-1,2]-L[m-2,1]
        # for Chebyshev T, use
        # L[m,1] = L[m-1,2]*sqrt(2)-L[m-2,1]
        # L[m,2] = L[m-1,3]+L[m-1,1]*sqrt(2)-L[m-2,2]
        # and go from j=3:m-1
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


function jacobioperator(a,b,t0,t1)
    Δ=ToeplitzOperator([t1],[t0,t1])
    B=jacobimatrix(a-t0,b-t1,length(a))
    Δ+CompactOperator(B)
end



immutable ToeplitzGivens <: BandedOperator{Float64}
    c::Float64
    s::Float64
end


function givenstail(t0::Real,t1::Real)
    @assert t0^2-4t1^2≥0
    s∞ = (t0 - sqrt(t0^2-4t1^2))/(2t1)
    c∞ = -sqrt(1-s∞^2)
    α = t1*c∞
    β = c∞*t0 - s∞*α
    l0 = (t0 + sqrt(t0^2-4t1^2))/2
    l1 = 2t1
    l2 = t1*s∞
    ToeplitzGivens(c∞,s∞),ToeplitzOperator([l1,l2],[l0]),α,β
end

#bandinds(T::ToeplitzGivens)=-ceil(Int,(-36-2log(abs(c)))/log(abs(s))),1
bandinds(T::ToeplitzGivens)=floor(Int,36/log(abs(T.s))),1

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

addentries!(T::ToeplitzGivens,A,kr::Range)=addentries!(ToeplitzOperator(T),A,kr)


function partialgivens(TG::ToeplitzGivens,m)
    T=ToeplitzOperator(TG)
    K=zeros(Float64,m-bandinds(T,1),m)
    neg,nonneg=T.negative,T.nonnegative
    for j=1:m-1
        if j > 1
            K[j-1,j]-=nonneg[2]
        end
        K[j,j]+=1-nonneg[1]
        for k=1:length(neg)
            K[k+j,j]-=neg[k]
        end
    end

    K[m-1,m]-=nonneg[2]
    c,s=TG.c,TG.s
    ret=c
    K[m,m]=c-nonneg[1]
    for k=1:length(neg)
        ret*=-s
        K[k+m,m]=ret-neg[k]
    end
    T+CompactOperator(K)
end


function ql(a,b,t0,t1)
    @assert t0^2>=4t1^2
    # The Givens rotations coming from infinity (with parameters c∞ and s∞) leave us with the almost triangular
    # a[n-1]  b[n-1]   0    0    0
    # b[n-1]   a[n]   t1    0    0
    #   0       α      β    0    0
    #   0      l2     l1   l0    0
    #   0       0     l2   l1   l0

    TQ,TL,α,β=givenstail(t0,t1)

    # Here we construct this matrix as L
    n = length(a)
    L = jacobimatrix(a,b,n+1)
    L[n,n+1] = t1
    #    L[n+1,n+2] = 0
    L[n+1,n+1]=β
    L[n+1,n]=α

    Q,L=tridql!(L)

    for k=1:size(Q,1)
        Q[k,k]-=1
    end
    for j=1:n+1
        L[j,j]-=TL.nonnegative[1]
        if j ≤ n
            L[j+1,j]-=TL.negative[1]
            if j ≤ n-1
                L[j+2,j]-=TL.negative[2]
            end
        end
    end

    partialgivens(TQ,n+1)*(I+CompactOperator(Q)),TL+CompactOperator(L)
end

t0=4.;t1=1.0; TQ,TL,α,β=givenstail(t0,t1)
Q,L=ql([1.,2.,3.],[5.,2.,1.],t0,t1)
J=jacobioperator([1.,2.,3.],[5.,2.,1.],t0,t1)

t1=0.5; a,b=[-2.,0.],[0.5];
J=jacobioperator(a,b,t0,t1)
for k=1:3
    μ=eigvals([a[1] b[1]; b[1] a[2]])[1]
    Q,L=ql(a-μ,b,-μ,t1)
    LQ=L*Q
    a,b=Float64[LQ[k,k] for k=1:length(a)+1]+μ,  Float64[LQ[k,k+1] for k=1:length(a)]
end
a,b



x=Fun((x,y)->x,Chebyshev()⊗Chebyshev())
y=Fun((x,y)->y,Chebyshev()⊗Chebyshev())


J=jacobioperator([0.,0.],[0.5],0.,0.5)
Jx=J⊗I
Jy=I⊗J

bandinds(Jx)
Jx[1:10,1:10][10,10]


convert(Matrix{Float64},(Jx^2+Jy^2)[1:10,1:10])

convert(Matrix{Float64},Jx[1:10,1:10])|>eigvals

convert(Matrix{Float64},Jy[1:10,1:10])|>eigvals

convert(Matrix{Float64},Jx[1:5,1:5])

Jx[1:5,1:5]


Jx=Multiplication(x,space(x))


Jx[1:10,1:10]


jacobioperator(a,b,0.,t1)[1:10,1:10]



μ
J2[1:10,1:10]|>chopm

J[1:10,1:10]|>full|>eigvals

LQ[1:10,1:10]|>chopm

J2[1:10,1:10]|>chopm

J=jacobioperator(a,b,t0,t1)


LQ=L*Q




TL[1:10,1:10]
using SO
(Q*L-J)[1:10,1:10]|>chopm

L*Q|>ApproxFun.introspect

J[1:10,1:10]

(L*Q)[5,5]

J[1:10,1:10]

(Q*L-J)[1:10,1:10]|>chopm

L[1:10,1:10]


full(TL[1:10,1:10])

full(TQ[1:10,1:10])*full(TL[1:10,1:10])






Q[1:n,1:n]*L[1:n,1:n]|>chopm

J=jacobioperator([1.,2.,3.],[5.,2.,1.],2.,0.5)

(Q'Q-I)[1:10,1:10]|>full|>norm

full(Q[1:10,1:10])'*full(J[1:10,1:10])

L.ops[1][1:10,1:10]
(Q'*J)[1:10,1:10]|>chopm

TL[500,500]/2

J[1:10,1:10]

QM*QM'-eye(n)|>chopm
QM=Q[1:n,1:n]|>full




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






## T vs U


function Lmatrix2(a,b,N)
    n = length(a)
    @assert n-length(b)==1
    bext = [b; .5]
    L = zeros(N,N)
    L[1,1] = 1
    L[2,1] = -a[1]/bext[1]

    L[2,2] = 1/(sqrt(2)*bext[1])

    # the generic case.
    for i = 3:n+1
        L[i,1] = (L[i-1,2]/sqrt(2)-a[i-1]*L[i-1,1]-bext[i-2]*L[i-2,1])/bext[i-1]
        L[i,2] = (L[i-1,3]/2+L[i-1,1]/sqrt(2)-a[i-1]*L[i-1,2]-bext[i-2]*L[i-2,2])/bext[i-1]

        for j = 3:i
            L[i,j] = (L[i-1,j+1]/2+L[i-1,j-1]/2-a[i-1]*L[i-1,j]-bext[i-2]*L[i-2,j])/bext[i-1]
        end
    end
    # this case is where b[m],b[m-1],b[m-2] = 1/2, and a[m],a[m-1] = 0, like Chebyshev
    for m = n+2:N
        L[m,1] = L[m-1,2]*sqrt(2)-L[m-2,1]
        L[m,2] = L[m-1,3]+L[m-1,1]*sqrt(2)-L[m-2,2]
        for j = 3:m-1
            L[m,j] = L[m-1,j+1]-L[m-2,j]+L[m-1,j-1]
        end
        L[m,m] = -L[m-2,m] + L[m-1,m-1]
    end
    L
end


1/sqrt(2.)

Conversion(C

Lmatrix([0.,0.],[0.2],20)

Lmatrix2([0.,0.],[0.5],20) |>inv
Lmatrix([0.,0.],[.5],20)

Lmatrix3([0.,0.],[1.],20)



function Lmatrix3(a,b,N)
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
