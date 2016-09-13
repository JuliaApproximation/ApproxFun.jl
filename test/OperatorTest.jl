using ApproxFun, Base.Test, Compat
    import ApproxFun:Multiplication,InterlaceOperator

import Compat.view


# test fast copy is consistent with getindex


C=ToeplitzOperator([1.,2.,3.],[4.,5.,6.])

@test_approx_eq full(C[1:5,1:5])    [4.0 5.0 6.0 0.0 0.0
                                     1.0 4.0 5.0 6.0 0.0
                                     2.0 1.0 4.0 5.0 6.0
                                     3.0 2.0 1.0 4.0 5.0
                                     0.0 3.0 2.0 1.0 4.0]

C=Conversion(Ultraspherical(1),Ultraspherical(2))

@test_approx_eq full(C[1:5,1:5])     [1.0 0.0 -0.3333333333333333 0.0  0.0
                                      0.0 0.5  0.0               -0.25 0.0
                                      0.0 0.0  0.3333333333333333 0.0 -0.2
                                      0.0 0.0  0.0                0.25 0.0
                                      0.0 0.0  0.0                0.0  0.2]




for M in (ToeplitzOperator([1.,2.,3.],[4.,5.,6.]),
                HankelOperator([1.,2.,3.,4.,5.,6.,7.]),
                Conversion(Ultraspherical(1),Ultraspherical(2)),
                Multiplication(Fun([1.,2.,3.],Chebyshev()),Chebyshev()))
    B=M[1:5,1:5]

    for k=1:5,j=1:5
        @test_approx_eq B[k,j] M[k,j]
    end

    @test_approx_eq M[1:5,1:5][2:5,1:5] M[2:5,1:5]
    @test_approx_eq M[1:5,2:5] M[1:5,1:5][:,2:end]
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

Q=integrate(d)

@test norm((Q+I)*f-(integrate(f)+f)) < 100eps()
@test norm((Q)*f-(integrate(f))) < 100eps()

x=Fun(identity)
X=Multiplication(x,space(x))
d=Interval()


A=Conversion(Chebyshev(d),Ultraspherical(2,d))


@test norm(A\Fun(x.*f,rangespace(A))-(x.*f)) < 100eps()

@test norm((Conversion(Chebyshev(d),Ultraspherical(2,d))\(D^2*f))-f'') < 100eps()

@test norm(X*f-(x.*f)) < 100eps()

A=Conversion(Chebyshev(d),Ultraspherical(2,d))*X
@test norm((A*f.coefficients).coefficients-coefficients(x.*f,rangespace(A))) < 100eps()


## Special functions

x=Fun(identity)
@test norm(cos(x)-Fun(cos))<10eps()
@test norm(sin(x)-Fun(sin))<10eps()
@test norm(exp(x)-Fun(exp))<10eps()
@test norm(sin(x)./x-Fun(x->sinc(x/π)))<100eps()


P=ApproxFun.PermutationOperator([2,1])
@test_approx_eq P[1:4,1:4] [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]



## Periodic


d=PeriodicInterval(0.,2π)
a=Fun(t-> 1+sin(cos(10t)),d)
D=Derivative(d)
L=D+a
f=Fun(t->exp(sin(t)),d)
u=L\f

@test norm(L*u-f) < 100eps()

d=PeriodicInterval(0.,2π)
a1=Fun(t->sin(cos(t/2)^2),d)
a0=Fun(t->cos(12sin(t)),d)
D=Derivative(d)
L=D^2+a1*D+a0

f=Fun([1,2,3,4,5],space(a1))

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

f=Fun(exp)
@test_approx_eq (A.ops[end]*f)(0.1) f'(0.1)
@test_approx_eq ((x*D)*f)(0.1) 0.1*f'(0.1)
@test_approx_eq (A*f)(0.1) f'(0.1)+0.1*f''(0.1)
@test_approx_eq (B*f)(0.1) f'(0.1)+0.1*f''(0.1)
@test_approx_eq (C*f)(0.1) f'(0.1)+0.1*f''(0.1)

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
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]


S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos);lneumann(S)])
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]



S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())])
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]

S=Chebyshev()
io=ApproxFun.InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos)])
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]


S=Chebyshev()
io=ApproxFun.InterlaceOperator([Derivative(Chebyshev());InterlaceOperator(dirichlet(S))])
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]


S=Chebyshev()
io=ApproxFun.InterlaceOperator([Derivative(Chebyshev())+Fun(cos);InterlaceOperator(dirichlet(S))])
co=cache(io)
@test co[1:100,1:100] == io[1:100,1:100]
@test co[1:100,1:100] == io[1:100,1:100]
@test co[200:300,200:300] == io[1:300,1:300][200:300,200:300]



## Reverse


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

let C=cache(A)
    ApproxFun.resizedata!(C,5,:)
    ApproxFun.resizedata!(C,10,:)
    @test norm(C.data-A[1:10,1:39]) == 0
end
