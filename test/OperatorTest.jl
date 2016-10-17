using ApproxFun, Base.Test, Compat
    import ApproxFun:Multiplication,InterlaceOperator

import Compat.view


# test row/colstarts

function functionaltest(A)
    @test ApproxFun.rowstart(A,1) == 1
    @test ApproxFun.colstop(A,1) == 1
    @test A[1:10]' == A[1,1:10]
    @test A[1:10][3:10] == A[3:10]
    @test A[1:10] == [A[k] for k=1:10]

    co=cache(A)
    @test co[1:100] == A[1:100]
    @test co[1:100] == A[1:100]
    @test co[200:300] == A[1:300][200:300] == A[200:300]
end

function infoperatortest(A)
    B=A[1:5,1:5]

    for k=1:5,j=1:5
        @test_approx_eq B[k,j] A[k,j]
    end

    @test_approx_eq A[1:5,1:5][2:5,1:5] A[2:5,1:5]
    @test_approx_eq A[1:5,2:5] A[1:5,1:5][:,2:end]

    @test isfinite(ApproxFun.colstart(A,1)) && ApproxFun.colstart(A,1) > 0
    @test isfinite(ApproxFun.rowstart(A,1)) && ApproxFun.colstart(A,2) > 0

    co=cache(A)
    @test co[1:100,1:100] == A[1:100,1:100]
    @test co[1:100,1:100] == A[1:100,1:100]
    @test co[200:300,200:300] == A[1:300,1:300][200:300,200:300] == A[200:300,200:300]

    let C=cache(A)
        ApproxFun.resizedata!(C,5,:)
        ApproxFun.resizedata!(C,10,:)
        @test norm(C.data[1:10,1:C.datasize[2]]-A[1:10,1:C.datasize[2]]) ≤ eps()
    end
end


function almostbandedoperatortest(A)
    infoperatortest(A)
    @test isfinite(ApproxFun.colstop(A,1))
end

function bandedoperatortest(A)
    almostbandedoperatortest(A)
    @test isfinite(ApproxFun.rowstop(A,1))
end


functionaltest(Evaluation(Chebyshev(),0.1,1)-Evaluation(Chebyshev(),0.1,1))

# test fast copy is consistent with getindex


C=ToeplitzOperator([1.,2.,3.],[4.,5.,6.])

bandedoperatortest(C)

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




for M in (HankelOperator([1.,2.,3.,4.,5.,6.,7.]),
            Multiplication(Fun([1.,2.,3.],Chebyshev()),Chebyshev()))
    bandedoperatortest(M)
end



d=Interval(-10.,5.);
S=Chebyshev(d)


@test norm(Fun(Fun(Fun(exp,S),Ultraspherical(2,d)),S)-Fun(exp,S)) < 100eps()


@test_approx_eq copy(@compat view(Derivative(Ultraspherical(1)),1:2,1:2))[1,2] Derivative(Ultraspherical(1))[1,2]
@test_approx_eq exp(0.1) (Derivative()*Fun(exp,Ultraspherical(1)))(0.1)


T=Float64
f=Fun(exp)
d=domain(f)
D=Derivative(d)

Q=Integral(d)

bandedoperatortest(Q)

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

bandedoperatortest(A)

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

bandedoperatortest(D)
bandedoperatortest(L)

f=Fun(t->exp(sin(t)),d)
u=L\f

@test norm(L*u-f) < 100eps()

d=PeriodicInterval(0.,2π)
a1=Fun(t->sin(cos(t/2)^2),d)
a0=Fun(t->cos(12sin(t)),d)
D=Derivative(d)
L=D^2+a1*D+a0

bandedoperatortest(L)

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
bandedoperatortest(x*D)

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
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev());lneumann(S)])

almostbandedoperatortest(io)


S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos);lneumann(S)])

almostbandedoperatortest(io)



S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())])
almostbandedoperatortest(io)

S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos)])
almostbandedoperatortest(io)


S=Chebyshev()
io=ApproxFun.InterlaceOperator([Derivative(Chebyshev());InterlaceOperator(dirichlet(S))])
almostbandedoperatortest(io)

S=Chebyshev()
io=ApproxFun.InterlaceOperator([Derivative(Chebyshev())+Fun(cos);InterlaceOperator(dirichlet(S))])
almostbandedoperatortest(io)



## Reverse


bandedoperatortest(ApproxFun.Reverse(Chebyshev()))
bandedoperatortest(ApproxFun.ReverseOrientation(Chebyshev()))

@test ApproxFun.Reverse(Chebyshev())*Fun(exp) ≈ Fun(x->exp(-x))
@test ApproxFun.ReverseOrientation(Chebyshev())*Fun(exp) ≈ Fun(exp,[1,-1])


@test norm(ApproxFun.Reverse(Fourier())*Fun(t->cos(cos(t-0.2)-0.1),Fourier()) - Fun(t->cos(cos(-t-0.2)-0.1),Fourier())) < 10eps()
@test norm(ApproxFun.ReverseOrientation(Fourier())*Fun(t->cos(cos(t-0.2)-0.1),Fourier()) - Fun(t->cos(cos(t-0.2)-0.1),Fourier([π,-π]))) < 10eps()




## Newton iteration bug
S=Chebyshev([0.,7.])

ω=2π


N = u->Any[Fun(u(0.)-0.1);Fun(u(ω)-u(0.));Fun(u'(ω)-u'(0.));u''+u+u^3]

u=0.1Fun(cos,S)

D=Derivative(S)

Z=ApproxFun.ZeroOperator(ApproxFun.ConstantSpace())

A=ApproxFun.interlace([Z                      Evaluation(S,0);
                     u'(ω)    Evaluation(S,ω)-Evaluation(S,0);
                     u''(ω)   Evaluation(S,ω,1)-Evaluation(S,0,1);
                      0         D^2+I+3u^2])

almostbandedoperatortest(A)
