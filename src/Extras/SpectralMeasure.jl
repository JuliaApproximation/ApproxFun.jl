using ApproxFun
    import ApproxFun:BandedOperator,ToeplitzOperator,tridql!,bandinds,DiracSpace, plot

joukowsky(z)=.5*(z+1./z)

function spectralmeasure(a,b)
    T,K=tkoperators(a,b)
    Tfun = Fun([T[1,1];T.negative],Taylor)
    eigs=sort!(real(map!(joukowsky,filter!(z->abs(z)<1,complexroots(Tfun)))))

    if isempty(eigs)
        L=T+K
        coeffs = L\[1]
        Fun(coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))
    else
        Q=Array(BandedOperator{Float64},0)

        for k=1:length(eigs)
            t1,t0=0.5,eigs[k]

            Q1,L1=ql(a-t0,b,-t0,t1)
            push!(Q,Q1)

            LQ=L1*Q1
            a=Float64[LQ[k,k] for k=2:length(a)+1]+t0;
            b=Float64[LQ[k,k+1] for k=2:length(a)];
            eigs[k]=LQ[1,1]+t0
        end

        L2=contspectraltransform(a,b)

        q0=[1.0]
        m=length(Q) #number of discrete
        for k=1:m
            q0=[q0[1:k-1];Q[k]'*q0[k:end]]
        end
        coeffs=L2\q0[m+1:end]
        c=q0[m+1]
        Fun([q0[1:m].^2;(2c/π)*coeffs],DiracSpace(JacobiWeight(0.5,0.5,Ultraspherical{1}()),eigs))
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

function tkoperators(a,b)
    @assert length(a)-length(b)==1
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

    # the generic case.
    for i = 3:n+1
        L[i,1] = (L[i-1,2]/2-a[i-1]*L[i-1,1]-bext[i-2]*L[i-2,1])/bext[i-1]
        for j = 2:i
            L[i,j] = (L[i-1,j+1]/2+L[i-1,j-1]/2-a[i-1]*L[i-1,j]-bext[i-2]*L[i-2,j])/bext[i-1]
        end
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
    l0 = (t0 + sqrt(t0^2-4t1^2))/2
    if s∞^2 > 1
        s∞ = (t0 + sqrt(t0^2-4t1^2))/(2t1)
        l0 = (t0 - sqrt(t0^2-4t1^2))/2
    end
    c∞ = -sqrt(1-s∞^2)
    α = t1*c∞
    β = c∞*t0 - s∞*α
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

sum(μ)

μ.coefficients[1:4]

n=20;a=zeros(n);b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


n=8;a=1/([1:n;]*sqrt(n));b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)



n=20;a=zeros(n);b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


n=11;a=0.2rand(n);b=rand(n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)




n=14;a=(rand(n)-0.5);b=fill(0.5,n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)



n=20;a=0.5cos([1:n]);b=fill(0.5,n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


α=44;β=44;a=[α/50+0.02145238,0.];b=[β/20+0.012412];
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)

values(μ)

ApproxFun.plot(μ)
μ
T,K = tkoperators(a,b)
T[1:10,1:10]
K[1:10,1:10]

