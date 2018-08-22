using ApproxFun, Test, StaticArrays
    import ApproxFun: testbandedbelowoperator, testbandedoperator, testspace, testtransforms, Vec,
                        maxspace, NoSpace, hasconversion, testfunctional

@testset "Jacobi" begin
    @test ApproxFun.jacobip(0:5,2,0.5,0.1) ≈ [1.,0.975,-0.28031249999999996,-0.8636328125,-0.0022111816406250743,0.7397117980957031]

    testspace(Jacobi(.5,2.);haslineintegral=false)

    f=Fun(exp,Jacobi(.5,2.))
    @test f(.1) ≈ exp(.1)

    f=Fun(x->cos(100x),Jacobi(.5,2.124),500)
    @test f(.1) ≈ cos(100*.1)


    sp=Jacobi(.5,2.124)
    @time f=Fun(exp,sp)
    sp2=Jacobi(1.5,2.124)
    f2=Fun(exp,sp2)
    sp3=Jacobi(1.5,3.124)
    f3=Fun(exp,sp3)
    sp4=Jacobi(2.5,4.124)
    f4=Fun(exp,sp4)
    @test norm((Fun(f,sp2)-f2).coefficients)<100eps()
    @test norm((Fun(f,sp3)-f3).coefficients)<100eps()
    @test norm((Fun(f,sp4)-f4).coefficients)<200eps()





    m=20
    @time testtransforms(JacobiWeight(0.,m,Jacobi(0.,2m+1)))
    f=Fun(x->((1-x)/2).^m.*exp(x),JacobiWeight(0.,m,Jacobi(0.,2m+1)))
    @test abs(f(.1)-(x->((1-x)/2).^m.*exp(x))(.1))<10eps()


    m=10
    @time f=Fun(x->besselj(m,m*(1-x)),JacobiWeight(0.,m,Jacobi(0.,2m+1)))
    @test f(0.) ≈ besselj(m,m)






    ## Conversion

    testtransforms(Jacobi(-0.5,-0.5))
    @test norm(Fun(Fun(exp),Jacobi(-.5,-.5))-Fun(exp,Jacobi(-.5,-.5))) < 300eps()

    x=Fun(identity)
    ri=0.5/(1-x)
    @test ((1-x)/2 .* Fun(exp,JacobiWeight(0.,0.,Jacobi(0.,1.))))(.1) ≈ (1-.1)./2*exp(.1)


    @test ((1-x)/2 .* Fun(exp,JacobiWeight(0.,0.,Jacobi(0.,1.))))(.1) ≈ (1-.1)./2*exp(.1)


    @test (ri*Fun(exp,JacobiWeight(0.,0.,Jacobi(0.,1.))))(.1) ≈ .5/(1-.1)*exp(.1)


    ## Derivative

    D=Derivative(Jacobi(0.,1.,Segment(1.,0.)))
    @time testbandedoperator(D)

    S=JacobiWeight(0.,0.,Jacobi(0.,1.,Segment(1.,0.)))
    D=Derivative(S)
    testbandedoperator(D)

    f=Fun(exp,domainspace(D))
    @test (D*f-f).coefficients|>norm < eps(100000.)
    @test (f'-f).coefficients|>norm < eps(100000.)
    @test (D^2*f-f).coefficients|>norm < eps(100000000.)
    @test (D*(D*f)-f).coefficients|>norm < eps(100000000.)





    ### Jacobi multiplication

    x=Fun(identity,Jacobi(0.,0.))
    f=Fun(exp,Jacobi(0.,0.))

    @test (x*f)(.1) ≈ .1exp(.1)

    x=Fun(identity,Jacobi(12.324,0.123))
    f=Fun(exp,Jacobi(0.,0.))

    @test (x*f)(.1) ≈ .1exp(.1)


    x=Fun(identity,Jacobi(12.324,0.123))
    f=Fun(exp,Jacobi(0.590,0.213))

    @test (x*f)(.1) ≈ .1exp(.1)

    g=Fun(cos,Jacobi(12.324,0.123))
    f=Fun(exp,Jacobi(0.590,0.213))

    @test (g*f)(.1) ≈ cos(.1)*exp(.1)


    ## Jacobi integrate and sum

    testtransforms(Legendre(0..2))
    @test sum(Fun(exp,Legendre(0..2))) ≈ sum(Fun(exp,0..2))

    a=Arc(0.,.1,0.,π/2)
    g=Fun(exp,Legendre(a))

    @test sum(g) ≈ sum(Fun(exp,a))



    ## Test special derivative

    x=Fun()
    f=exp(x)*sqrt(1-x^2)
    D=Derivative(WeightedJacobi(.5,.5))

    testtransforms(WeightedJacobi(.5,.5))
    testbandedoperator(D)

    @time g=(D*Fun(f,domainspace(D)))
    @test f'(0.1) ≈ g(0.1)

    ## Test implementation of conversion between Chebyshev and Jacobi spaces using FastTransforms

    f = Fun(x->cospi(1000x))
    g = Fun(f,Legendre())
    h = Fun(g,Chebyshev())
    @test norm(f.coefficients-h.coefficients,Inf) < 100eps()
    @time h = Fun(h,Legendre())
    @test norm(g.coefficients-h.coefficients,Inf) < 1000eps()





    ## Derivative

    S=JacobiWeight(1,1,Ultraspherical(1))

    f=Fun(S,[1.,2.,3.])
    @test (Derivative(S,2)*f)(0.1) ≈ f''(0.1)


    ## == tests


    @test WeightedJacobi(0.1,0.2) == WeightedJacobi(0.1+eps(),0.2)

    # this tests a subspace bug
    f=Fun(WeightedJacobi(0.1,0.2),rand(10))  # convert to Legendre expansion

    g=(f|(2:ApproxFun.∞))

    @test ApproxFun.coefficients(g.coefficients,space(g),ApproxFun.canonicalspace(g))[1] ==0.
    @test norm((Fun(g,space(f))|(2:ApproxFun.∞)-g).coefficients) < 10eps()


    ## Check conversion for non-compatible paramters
    S=Jacobi(1.2,0.1)
    x=Fun()


    p=(S,k)->Fun(S,[zeros(k);1.])
    n=1;
    @test norm(x*p(S,n-1)-(ApproxFun.recα(Float64,S,n)*p(S,n-1) + ApproxFun.recβ(Float64,S,n)*p(S,n))) < 10eps()



    # Log with squareroot singularities

    a = 1.0; b = 2.0+im
    d = Segment(a,b)
    z = Fun(d)

    f = real(exp(z) / (sqrt(z-a)*sqrt(b-z)))
    S=space(f)
    x=4.0+2im;
    @test linesum(f*log(abs(x-z))) ≈ 13.740676344264614



    # Line sum for legendre

    x=Fun(Legendre())
    @test sum(x+1) ≈ linesum(x+1)

    x=Fun(Legendre(2..1))
    @test sum(x+1) ≈ -linesum(x+1)

    x=Fun(1..1+im)
    @test sum(x+1) ≈ im*linesum(x+1)

    x=Fun(Legendre(Segment(1,1+im)))
    @test sum(x+1) ≈ im*linesum(x+1)

    x=Fun(Legendre(Segment(im,1)))
    @test sum(x+1) ≈ (1-im)/sqrt(2)*linesum(x+1)


    ## Laguerre and Hermite

    f=Fun(Laguerre(0.),[1,2,3])
    @test f(0.1) ≈ 5.215


    f = Fun(Laguerre(0.1),fill(1.0,100))
    @test f(0.2) ≈ 8.840040924281498


    @test (Derivative(Laguerre(0.1))*f)(0.2) ≈ -71.44556705957386
    f = Fun(Laguerre(0.2),fill(1.0,100))
    @test (Derivative(Laguerre(0.2))*f)(0.3) ≈ -137.05785783078218


    @test (Conversion(Laguerre(0.2),Laguerre(1.2))*f)(0.1) ≈ f(0.1)
    @test (Conversion(Laguerre(0.2),Laguerre(2.2))*f)(0.1) ≈ f(0.1)



    f=Fun(LaguerreWeight(0.,Laguerre(0.1)),fill(1.0,100))
    @test f'(0.2) ≈ -65.7322962859456


    B=Evaluation(LaguerreWeight(0.,Laguerre(0.1)),false)
    @test B*f ≈ 151.53223385808576


    x=Fun(Laguerre(0.0))
    S=WeightedLaguerre(0.0)
    D=Derivative(S)
    u=[ldirichlet();D^2-x]\[airyai(0.0);0.0]
    @test u(1.0) ≈ airyai(1.0)

    w = Fun(WeightedLaguerre(), [1.0])
    @test sum(w) == 1
    t = Fun(identity, space(w))
    @test t(10.0) == 10.0
    M = Multiplication(t^2+1, space(w))
    (M \ w)(1.0) ≈ exp(-1)/2
    #

    ## Test vector valued case

    f=Fun((x,y)->real(exp(x+im*y)),Legendre(Vec(0.,0.)..Vec(1.,1.)))
    @test f(0.1,0.1) ≈ real(exp(0.1+0.1im))



    ## Test integer, float mixed

    C=Conversion(Legendre(),Jacobi(1,0))
    testbandedoperator(C)


    @testset "Addition of piecewise Legendre bug" begin
        f = Fun(exp,Legendre())
        f1 = Fun(exp,Legendre(-1..0))
        f2 = Fun(exp,Legendre(0..1))
        fp = f1+f2
        @test space(fp) isa PiecewiseSpace
        @test fp(0.1) ≈ exp(0.1)
        @test fp(0.) ≈ exp(0.)
        @test fp(-0.1) ≈ exp(-0.1)
    end



    @testset "JacobiWeight cumsum bug Issue #557" begin
        x = Fun(0.0..1.0)
        ν = 2
        @time f = x^(ν/2-1) * exp(-x/2) # 0.05s
        @test cumsum(f)' ≈ f
        @test cumsum(f)(1.0) ≈ 0.7869386805747332 # Mathematic

        x = Fun(0.0..Inf)
        ν = 2
        @time f = x^(ν/2-1) * exp(-x/2) # 0.05s
        @test cumsum(f)' ≈ f
        @test cumsum(f)(1.0) ≈ 0.7869386805747332
    end


    @testset "Jacobi–Chebyshev conversion tests" begin
        a,b = (Jacobi(-0.5,-0.5), Legendre())
        @test maxspace(a,b) == NoSpace()
        @test union(a,b) == a
        @test !hasconversion(a,b)

        a,b = (Chebyshev(), Legendre())
        @test maxspace(a,b) == NoSpace()
        @test union(a,b) == Jacobi(-0.5,-0.5)
        @test !hasconversion(a,b)
    end

    @testset "Definite integral tests" begin
        B = DefiniteIntegral(WeightedJacobi(0,0))
        testfunctional(B)
        @test ApproxFun.rowstop(B,1) == 1
    end

    @testset "WeightedLaguerre cumsum" begin
        α = 2.7
        f = Fun(WeightedLaguerre(α), [1.0]);
        f = Fun(f, JacobiWeight(α,0,Chebyshev(0.0 .. Inf)));
        g = integrate(f)
        g(3.0) - cumsum(Fun(x -> f(x), 0..6))(3.0)
    end

    @testset "Definite integral tests" begin
        for S in (WeightedJacobi(0,0), JacobiWeight(0,0, Legendre(1.1..2.3)), Legendre())
            B = DefiniteIntegral(S)
            testfunctional(B)
            @test ApproxFun.rowstop(B,1) == 1
            B[1] == arclength(domain(S))
            f = Fun(exp, S)
            B*f == sum(Fun(exp,domain(S)))
        end
    end
end
