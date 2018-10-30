using ApproxFun, LinearAlgebra, SpecialFunctions, FastTransforms, Test
    import ApproxFun: testbandedoperator, testraggedbelowoperator, InterlaceOperator, testspace,
                        testbandedbelowoperator, testbandedblockbandedoperator, testfunctional, factor
## This includes extra tests that are too time consuming for Travis


include("runtests.jl")

@time include("FractionalTest.jl")

@testset "Full Operator" begin
    @time for M in (Multiplication(Fun(CosSpace(),[1.]),CosSpace()),
                    Multiplication(Fun(CosSpace(),[1.]),SinSpace()),
                    Multiplication(Fun(SinSpace(),[1.]),SinSpace()),
                    Multiplication(Fun(SinSpace(),[1.]),CosSpace()),
                    Derivative(SinSpace()),Derivative(CosSpace()))
          testbandedoperator(M)
    end

    S=Chebyshev()
    @time for io in (
            [Dirichlet(S);Derivative(Chebyshev());lneumann(S)],
            [Dirichlet(S);Derivative(Chebyshev())+Fun(cos);lneumann(S)],
            [Dirichlet(S);Derivative(Chebyshev())],
            [Dirichlet(S);Derivative(Chebyshev())+Fun(cos)],
            [Derivative(Chebyshev());Dirichlet(S)],
            [Derivative(Chebyshev())+Fun(cos);Dirichlet(S)])
        testraggedbelowoperator(io)
    end

<<<<<<< HEAD
d=ChebyshevInterval()
D=Derivative(d)
A=D^2-I
@time κ=nullspace(A)
@test length(κ) == 2
=======
    ## Newton iteration bug
    S=Chebyshev(0..7)
>>>>>>> abff326fa184c4021c60a8af5d7be726eccfbe54

    ω=2π


    N = u->Any[Fun(u(0.)-0.1);Fun(u(ω)-u(0.));Fun(u'(ω)-u'(0.));u''+u+u^3]

    u=0.1Fun(cos,S)

    D=Derivative(S)

    Z=ApproxFun.ZeroOperator(ApproxFun.ConstantSpace())


    testfunctional(Evaluation(S,0))
    testfunctional(Evaluation(S,ω)-Evaluation(S,0))
    testfunctional(Evaluation(S,ω,1)-Evaluation(S,0,1))

    A=[Z                      Evaluation(S,0);
                         u'(ω)    Evaluation(S,ω)-Evaluation(S,0);
                         u''(ω)   Evaluation(S,ω,1)-Evaluation(S,0,1);
                          0         D^2+I+3u^2]

    ApproxFun.backend_testinfoperator(A)

    DS=WeightedJacobi(0.1+1,0.2+1)
    D=Derivative(DS)[2:end,:]
    @time ApproxFun.testbandedoperator(D)
end

@testset "Full ODE" begin
    ## Null space
    d=ChebyshevInterval()
    D=Derivative(d)
    A=D^2-I
    @time κ=nullspace(A)
    @test length(κ) == 2

    c=[κ(0.);κ'(0.)]\[exp(0.),exp(0.)]
    u=(κ*c)[1]

    @test u(1.0) ≈ ℯ


    d=Interval(-50.,5.)
    x=Fun(identity,d)
    D=Derivative(d)
    @time u=nullspace(D^2-x)
    c=[u(leftendpoint(d)); u(rightendpoint(d))]\[airyai(d.a),airyai(d.b)]
    @test norm((u*c)[1]-Fun(airyai,d))<10000eps()


    ## constant forcing
    d = Interval(0.,50.)
    D = Derivative(d)
    t = Fun(identity,d)

    F = D^2 +.5D + I

    A= [ 0    ldirichlet(d);
         0    lneumann(d);
         0    rdirichlet(d);
        -1    F; ]

    @time u,x=A\[1.,0.,2.,0.]

    @test norm(F*x-u)<1000eps()


    @testset "Bessel" begin
        @time for ν in (1.,0.5,2.,3.5)
            println("        ν = $ν")
            S=JacobiWeight(-ν,0.,Chebyshev(0..1))
            D=Derivative(S)
            x=Fun(identity,domain(S))
            L=(x^2)*D^2+x*D+(x^2-ν^2);
            u=\([rdirichlet(S);rneumann(S);L],[bessely(ν,1.),.5*(bessely(ν-1.,1.)-bessely(ν+1.,1.)),0];
                        tolerance=1E-10)
            @test ≈(u(.1),bessely(ν,.1);atol=eps(1000000.)*max(abs(u(.1)),1))
            u=Fun(x->bessely(ν,x),S)
            @test ≈(u(.1),bessely(ν,.1);atol=eps(10000.)*max(abs(u(.1)),1))
            u=Fun(x->besselj(ν,x),S)
            @test ≈(u(.1),besselj(ν,.1);atol=eps(10000.)*max(abs(u(.1)),1))
        end

        @time for ν in (1.,0.5,0.123,3.5)
            println("        ν = $ν")
            S=JacobiWeight(ν,0.,Chebyshev(0..1))
            D=Derivative(S)
            x=Fun(identity,domain(S))
            L=(x^2)*D^2+x*D+(x^2-ν^2);

            u=\([rdirichlet(S);rneumann(S);L],[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.)),0];
                        tolerance=1E-10)
            @test ≈(u(.1),besselj(ν,.1);atol=eps(1000000.)*max(abs(u(.1)),1))
            u=Fun(x->besselj(ν,x),S)
            @test ≈(u(.1),besselj(ν,.1);atol=eps(10000.)*max(abs(u(.1)),1))
        end
    end
end

@testset "Full Jacobi" begin
    sp = Jacobi(.5,2.124)
    f = Fun(exp,sp)
    sp2 = Jacobi(1.5,2.124)
    M = Multiplication(f,sp2)
    @time testbandedoperator(M)


    ## Legendre conversions
    testspace(Ultraspherical(1); haslineintegral=false)
    testspace(Ultraspherical(2); haslineintegral=false)
    # minpoints is a tempory fix a bug
    @time testspace(Ultraspherical(1//2); haslineintegral=false, minpoints=2)
    @test norm(Fun(exp,Ultraspherical(1//2))-Fun(exp,Jacobi(0,0))) < 100eps()

    C=Conversion(Jacobi(0,0),Chebyshev())
    @time testbandedbelowoperator(C)
    @test norm(C*Fun(exp,Jacobi(0,0))  - Fun(exp)) < 100eps()


    C=Conversion(Ultraspherical(1//2),Chebyshev())
    @time testbandedbelowoperator(C)
    @test norm(C*Fun(exp,Ultraspherical(1//2))  - Fun(exp)) < 100eps()



    C=Conversion(Chebyshev(),Ultraspherical(1//2))
    @time testbandedbelowoperator(C)
    @test norm(C*Fun(exp)-Fun(exp,Legendre())) < 100eps()


    C=Conversion(Chebyshev(),Jacobi(0,0))
    @time testbandedbelowoperator(C)
    @test norm(C*Fun(exp)  - Fun(exp,Jacobi(0,0))) < 100eps()


    C=Conversion(Chebyshev(),Jacobi(1,1))
    @time testbandedbelowoperator(C)
    @test norm(C*Fun(exp) - Fun(exp,Jacobi(1,1))) < 100eps()


    C=Conversion(Ultraspherical(1//2),Ultraspherical(1))
    @time testbandedbelowoperator(C)

    λ1 = ApproxFun.order(domainspace(C))
    λ2 = ApproxFun.order(rangespace(C))

    # test against version that doesn't use lgamma
    Cex = Float64[(if j ≥ k && iseven(k-j)
            gamma(λ2)*(k-1+λ2)/(gamma(λ1)*gamma(λ1-λ2))*
                (gamma((j-k)/2+λ1-λ2)/gamma((j-k)/2+1))*
                (gamma((k+j-2)/2+λ1)/gamma((k+j-2)/2+λ2+1))
        else
            0.0
        end) for k=1:20,j=1:20]

    @test norm(Cex - C[1:20,1:20]) < 100eps()

    @test norm(C*Fun(exp,Ultraspherical(1//2))-Fun(exp,Ultraspherical(1))) < 100eps()

    C=Conversion(Jacobi(0,0),Ultraspherical(1))
    testbandedbelowoperator(C)
    @test norm(C*Fun(exp,Jacobi(0,0))-Fun(exp,Ultraspherical(1))) < 100eps()


    C=Conversion(Ultraspherical(1),Jacobi(0,0))
    testbandedbelowoperator(C)
    @test norm(C*Fun(exp,Ultraspherical(1))-Fun(exp,Jacobi(0,0))) < 100eps()
end

@testset "Full multivariate" begin
    ## ProductFun
    u0   = ProductFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-0.1)^2)+.5exp(-30(x+0.5)^2-40(y+0.2)^2))


    @test values(u0)-values(u0|>LowRankFun)|>norm < 1000eps()
    @test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

    ##TODO: need to do adaptive to get better accuracy
    @test sin(u0)(.1,.2)-sin(u0(.1,.2))|>abs < 10e-4


    F = LowRankFun((x,y)->hankelh1(0,10abs(y-x)),Chebyshev(1.0..2.0),Chebyshev(1.0im..2.0im))

    @test F(1.5,1.5im) ≈ hankelh1(0,10abs(1.5im-1.5))


    ## Periodic
    f=LowRankFun((x,y)->cos(x)*sin(y),PeriodicSegment(),PeriodicSegment())
    @test f(.1,.2) ≈ cos(.1)*sin(.2)

    f=LowRankFun((x,y)->cos(cos(x)+sin(y)),PeriodicSegment(),PeriodicSegment())
    @test f(.1,.2) ≈ cos(cos(.1)+sin(.2))
    @test norm(Float64[cos(cos(x)+sin(y)) for x=ApproxFun.vecpoints(f,1),y=ApproxFun.vecpoints(f,2)]-values(f))<10000eps()

    f=ProductFun((x,y)->cos(cos(x)+sin(y)),PeriodicSegment()^2)
    @test f(.1,.2) ≈ cos(cos(.1)+sin(.2))
    x,y=points(f)
    @test norm(Float64[cos(cos(x[k,j])+sin(y[k,j])) for k=1:size(f,1),j=1:size(f,2)]-values(f))<10000eps()

    d=PeriodicSegment()^2
    f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
    @test (transpose(f)-f|>coefficients|>norm)< 1000eps()

    ## Functional*Fun

<<<<<<< HEAD
d=ChebyshevInterval()
B=ldirichlet(d)
f=ProductFun((x,y)->cos(cos(x)*sin(y)),d^2)
=======
    d=Interval()
    B=ldirichlet(d)
    f=ProductFun((x,y)->cos(cos(x)*sin(y)),d^2)
>>>>>>> abff326fa184c4021c60a8af5d7be726eccfbe54

    @test norm(B*f-Fun(y->cos(cos(-1)*sin(y)),d))<20000eps()
    @test norm(f*B-Fun(x->cos(cos(x)*sin(-1)),d))<20000eps()

    ## matrix

<<<<<<< HEAD
f=Fun((x,y)->[exp(x*cos(y));cos(x*sin(y));2],ChebyshevInterval()^2)
@test f(0.1,0.2) ≈ [exp(0.1*cos(0.2));cos(0.1*sin(0.2));2]

f=Fun((x,y)->[exp(x*cos(y)) cos(x*sin(y)); 2 1],ChebyshevInterval()^2)
@test f(0.1,0.2) ≈ [exp(0.1*cos(0.2)) cos(0.1*sin(0.2));2 1]
=======
    f=Fun((x,y)->[exp(x*cos(y));cos(x*sin(y));2],Interval()^2)
    @test f(0.1,0.2) ≈ [exp(0.1*cos(0.2));cos(0.1*sin(0.2));2]

    f=Fun((x,y)->[exp(x*cos(y)) cos(x*sin(y)); 2 1],Interval()^2)
    @test f(0.1,0.2) ≈ [exp(0.1*cos(0.2)) cos(0.1*sin(0.2));2 1]
>>>>>>> abff326fa184c4021c60a8af5d7be726eccfbe54


    ## Cauchy fun

    f = Fun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),Line()^2)
    @test f(0.1,0.2) ≈ 1/(2π*(0.1^2+0.2^2+1)^(3/2))


    #TODO: improve tolerance
    f = LowRankFun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),JacobiWeight(2.,2.,Line())^2)
    @test ≈(f(0.1,0.2),1/(2π*(0.1^2+0.2^2+1)^(3/2));atol=1E-4)



    # 2d derivative (issue #346)
    @time let d = Chebyshev()^2
        f = Fun((x,y) -> sin(x) * cos(y), d)
        C=Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical(1)⊗Ultraspherical(1))
        @test (C*f)(0.1,0.2) ≈ f(0.1,0.2)
        Dx = Derivative(d, [1,0])
        f = Fun((x,y) -> sin(x) * cos(y), d)
        fx = Fun((x,y) -> cos(x) * cos(y), d)
        @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
        Dy = Derivative(d, [0,1])
        fy = Fun((x,y) -> -sin(x) * sin(y), d)
        @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
        L=Dx+Dy
        testbandedblockbandedoperator(L)

        @test (L*f)(0.2,0.3) ≈ (fx(0.2,0.3)+fy(0.2,0.3))

        B=ldirichlet(factor(d,1))⊗ldirichlet(factor(d,2))
        @test Number(B*f) ≈ f(-1.,-1.)

        B=Evaluation(factor(d,1),0.1)⊗ldirichlet(factor(d,2))
        @test Number(B*f) ≈ f(0.1,-1.)

        B=Evaluation(factor(d,1),0.1)⊗Evaluation(factor(d,2),0.3)
        @test Number(B*f) ≈ f(0.1,0.3)

        B=Evaluation(d,(0.1,0.3))
        @test Number(B*f) ≈ f(0.1,0.3)
    end
end

include("FullPDETest.jl")
println("Speed tests")
include("SpeedTest.jl")
include("SpeedODETest.jl")
include("SpeedPDETest.jl")

include("ReadmeTest.jl")
