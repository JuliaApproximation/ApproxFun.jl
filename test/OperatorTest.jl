using ApproxFun, Base.Test, Compat
    import ApproxFun:Multiplication,InterlaceOperator
    import Compat.view
    import ApproxFun: functionaltest, bandedoperatortest, raggedbelowoperatortest, infoperatortest


# test row/colstarts


functionaltest(Evaluation(Chebyshev(),0.1,1)-Evaluation(Chebyshev(),0.1,1))

# test fast copy is consistent with getindex


C=ToeplitzOperator([1.,2.,3.],[4.,5.,6.])

@time bandedoperatortest(C)

@test_approx_eq full(C[1:5,1:5])    [4.0 5.0 6.0 0.0 0.0
                                     1.0 4.0 5.0 6.0 0.0
                                     2.0 1.0 4.0 5.0 6.0
                                     3.0 2.0 1.0 4.0 5.0
                                     0.0 3.0 2.0 1.0 4.0]

C=Conversion(Ultraspherical(1),Ultraspherical(2))

bandedoperatortest(C)

@test_approx_eq full(C[1:5,1:5])     [1.0 0.0 -0.3333333333333333 0.0  0.0
                                      0.0 0.5  0.0               -0.25 0.0
                                      0.0 0.0  0.3333333333333333 0.0 -0.2
                                      0.0 0.0  0.0                0.25 0.0
                                      0.0 0.0  0.0                0.0  0.2]




@time for M in (HankelOperator([1.,2.,3.,4.,5.,6.,7.]),
            Multiplication(Fun([1.,2.,3.],Chebyshev()),Chebyshev()))
    bandedoperatortest(M)
end



d=Interval(-10.,5.);
S=Chebyshev(d)


@test norm(Fun(Fun(Fun(exp,S),Ultraspherical(2,d)),S)-Fun(exp,S)) < 100eps()


@test_approx_eq copy(@compat view(Derivative(Ultraspherical(1)),1:2,1:2))[1,2] Derivative(Ultraspherical(1))[1,2]
@test_approx_eq exp(0.1) (Derivative()*Fun(exp,Ultraspherical(1)))(0.1)


f=Fun(exp)
d=domain(f)
Q=Integral(d)
D=Derivative(d)

@time bandedoperatortest(Q)

@test norm((Q+I)*f-(integrate(f)+f)) < 100eps()
@test norm((Q)*f-(integrate(f))) < 100eps()

x=Fun(identity)
X=Multiplication(x,space(x))

bandedoperatortest(X)

d=Interval()


A=Conversion(Chebyshev(d),Ultraspherical(2,d))

bandedoperatortest(A)

@test norm(A\Fun(x.*f,rangespace(A))-(x.*f)) < 100eps()

@test norm((Conversion(Chebyshev(d),Ultraspherical(2,d))\(D^2*f))-f'') < 100eps()

@test norm(X*f-(x.*f)) < 100eps()

A=Conversion(Chebyshev(d),Ultraspherical(2,d))*X

@time bandedoperatortest(A)

@test norm((A*f.coefficients).coefficients-coefficients(x.*f,rangespace(A))) < 100eps()


## Special functions

x=Fun(identity)
@test norm(cos(x)-Fun(cos))<10eps()
@test norm(sin(x)-Fun(sin))<10eps()
@test norm(exp(x)-Fun(exp))<10eps()
@test norm(sin(x)./x-Fun(x->sinc(x/π)))<100eps()


P=ApproxFun.PermutationOperator([2,1])

bandedoperatortest(P)

@test_approx_eq P[1:4,1:4] [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]



## Periodic


d=PeriodicInterval(0.,2π)
a=Fun(t-> 1+sin(cos(10t)),d)
D=Derivative(d)
L=D+a

@time bandedoperatortest(D)
@time bandedoperatortest(Multiplication(a,Space(d)))


f=Fun(t->exp(sin(t)),d)
u=L\f

@test norm(L*u-f) < 100eps()

d=PeriodicInterval(0.,2π)
a1=Fun(t->sin(cos(t/2)^2),d)
a0=Fun(t->cos(12sin(t)),d)
D=Derivative(d)
L=D^2+a1*D+a0

@time bandedoperatortest(L)

f=Fun([1,2,3,4,5],space(a1))

bandedoperatortest(Multiplication(a0,Fourier([0.,2π])))

@test_approx_eq (Multiplication(a0,Fourier([0.,2π]))*f)(0.1)  (a0(0.1)*f(0.1))
@test_approx_eq ((Multiplication(a1,Fourier([0.,2π]))*D)*f)(0.1)  (a1(0.1)*f'(0.1))
@test_approx_eq (L.ops[1]*f)(0.1) f''(0.1)
@test_approx_eq (L.ops[2]*f)(0.1) a1(0.1)*f'(0.1)
@test_approx_eq (L.ops[3]*f)(0.1) a0(0.1)*f(0.1)
@test_approx_eq (L*f)(0.1) f''(0.1)+a1(0.1)*f'(0.1)+a0(0.1)*f(0.1)

f=Fun(t->exp(cos(2t)),d)
u=L\f

@test norm(L*u-f) < 1000eps()




## Check mixed

d=Interval()
D=Derivative(d)
x=Fun(identity,d)
A=D*(x*D)
B=D+x*D^2
C=x*D^2+D

bandedoperatortest(A)
bandedoperatortest(B)
bandedoperatortest(C)
@time bandedoperatortest(x*D)

f=Fun(exp)
@test_approx_eq (A.ops[end]*f)(0.1) f'(0.1)
@test_approx_eq ((x*D)*f)(0.1) 0.1*f'(0.1)
@test_approx_eq (A*f)(0.1) f'(0.1)+0.1*f''(0.1)
@test_approx_eq (B*f)(0.1) f'(0.1)+0.1*f''(0.1)
@test_approx_eq (C*f)(0.1) f'(0.1)+0.1*f''(0.1)

bandedoperatortest(A-B)
bandedoperatortest(B-A)
bandedoperatortest(A-C)

@test norm((A-B)[1:10,1:10]|>full)<eps()
@test norm((B-A)[1:10,1:10]|>full)<eps()
@test norm((A-C)[1:10,1:10]|>full)<eps()
@test norm((C-A)[1:10,1:10]|>full)<eps()
@test norm((C-B)[1:10,1:10]|>full)<eps()
@test norm((B-C)[1:10,1:10]|>full)<eps()



## Cached operator
@test cache(Derivative(Chebyshev(),2))[1,1]==0


S=Chebyshev()
D=Derivative(S)
@time for padding = [true,false]
  co=ApproxFun.CachedOperator(D,ApproxFun.RaggedMatrix(Float64[],Int[1],0),(0,0),domainspace(D),rangespace(D),bandinds(D),padding) #initialise with empty RaggedMatrix
  @test co[1:20,1:10] == D[1:20,1:10]
  @test size(co.data) == (20,10)
  ApproxFun.resizedata!(co,10,30)
  @test size(co.data)[2] == 30 && size(co.data)[1] ≥ 20
end

## Reverse


bandedoperatortest(ApproxFun.Reverse(Chebyshev()))
bandedoperatortest(ApproxFun.ReverseOrientation(Chebyshev()))

@test ApproxFun.Reverse(Chebyshev())*Fun(exp) ≈ Fun(x->exp(-x))
@test ApproxFun.ReverseOrientation(Chebyshev())*Fun(exp) ≈ Fun(exp,[1,-1])


@test norm(ApproxFun.Reverse(Fourier())*Fun(t->cos(cos(t-0.2)-0.1),Fourier()) - Fun(t->cos(cos(-t-0.2)-0.1),Fourier())) < 10eps()
@test norm(ApproxFun.ReverseOrientation(Fourier())*Fun(t->cos(cos(t-0.2)-0.1),Fourier()) - Fun(t->cos(cos(t-0.2)-0.1),Fourier([π,-π]))) < 10eps()





## Sub interval
f=Fun(exp)

D = Derivative(Chebyshev())
u = D[:,2:end] \ f
@test norm(u'-f) < 10eps()
@test_approx_eq u(0.1) exp(0.1)-f.coefficients[1]


u = D[1:end,2:end] \ f
@test_approx_eq u(0.1) exp(0.1)-f.coefficients[1]

u = D[1:ApproxFun.∞,2:ApproxFun.∞] \ f
@test_approx_eq u(0.1) exp(0.1)-f.coefficients[1]


DS=WeightedJacobi(0.1+1,0.2+1)
D=Derivative(DS)[2:end,:]

@time ApproxFun.bandedoperatortest(D)
