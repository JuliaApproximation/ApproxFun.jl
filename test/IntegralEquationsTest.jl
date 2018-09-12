using ApproxFun, IntervalSets, Random, Test
    import ApproxFun: testfunctional, testbandedbelowoperator, testbandedoperator


@testset "Integral equations" begin
    @time for S in (Fourier(Circle()),Laurent(Circle()),Taylor(Circle()),
                CosSpace(Circle()),JacobiWeight(-0.5,-0.5,Chebyshev()),
                JacobiWeight(-0.5,-0.5,Chebyshev(Segment(1.0,2.0+im))),
                JacobiWeight(0.5,0.5,Ultraspherical(1,Segment(1.0,2.0+im))))
        testfunctional(DefiniteLineIntegral(S))
    end

    # checks bug for
    dom = Segment(-1.0,1) ∪ Segment(1.0+im,1+2im) ∪ Segment(-2.0+im,-1+2im)


    ⨍ = DefiniteLineIntegral(dom)
    S = domainspace(⨍)
    @test ApproxFun.blockbandwidths(⨍) == (0,0)

    f = Fun(S,rand(20))

    @test DefiniteLineIntegral(component(dom,1))*component(f,1) ≈ linesum(component(f,1))
    @test DefiniteLineIntegral(component(dom,2))*component(f,2) ≈ linesum(component(f,2))
    @test DefiniteLineIntegral(component(dom,3))*component(f,3) ≈ linesum(component(f,3))

    @test ⨍*f ≈ linesum(f)

    #The first test checks the solution of the integral equation
    # u(x) + \int_{-1}^{+1} \frac{e^{y} u(y)}{\sqrt{1-y^2}} dy = f
    # on the interval -1..1.

    x = Fun(identity)
    w = 1/sqrt(1-x^2)
    d = domain(x)

    Σ = DefiniteIntegral(d)
    @test bandwidths(Σ) == (0,0)

    @test domainspace(Σ) ==
        JacobiWeight{Chebyshev{ChebyshevInterval{Float64}},ChebyshevInterval{Float64},Float64,Float64}(-0.5,-0.5,Chebyshev())

    L = I+Σ[exp(x)*w]
    testbandedoperator(L)
    usol = sin(2x)
    f = L*usol

    @time u = L\f
    @test norm(u-usol) <= 10eps()


    #The second test checks the solution of the integro-differential equation
    # u'(x) + x u(x) + \int_{-2}^{+2} sin(y-x) u(y) \sqrt{4-y^2} dy = f
    # on the interval -2..2, with u(-2) = 1.

    x = Fun(identity,-2..2)
    w = sqrt(4-x^2)
    d = domain(x)

    D = Derivative(d)
    B = ldirichlet(d)
    Σ = DefiniteIntegral(.5,.5,d)

    @test domainspace(Σ) ==
        JacobiWeight{Ultraspherical{Int,ClosedInterval{Int},Float64},ClosedInterval{Float64},Float64,Float64}(.5,.5,Ultraspherical(1,d))

    K = LowRankFun((x,y)->sin(y-x)*w(y),Ultraspherical(1,d),domainspace(Σ))

    L = D+x+Σ[K]
    usol = cospi(20x)
    f = L*usol
    @time u = [B;L]\[1.;f]


    @test norm(u-usol) ≤ 200eps()


    Σ = DefiniteIntegral()

    f1 = Fun(t->cos(cos(t)),-π..π)
    f = Fun(t->cos(cos(t)),Laurent(-π..π))

    @test sum(f1) ≈ Σ*f

    f1 = Fun(t->cos(cos(t))/t,Laurent(Circle()))
    f2 = Fun(t->cos(cos(t))/t,Fourier(Circle()))
    @test Σ*f1 ≈ Σ*f2

    f1=Fun(t->cos(cos(t)),Laurent(-π..π))
    f2=Fun(t->cos(cos(t)),Fourier(-π..π))
    @test Σ*f1 ≈ Σ*f2


    ## test over arcs


    d = exp(im*Segment(0.1,0.2))
    x = Fun(d)
    @time w=1/(sqrt(abs(first(d)-x))*sqrt(abs(last(d)-x)))

    @test linesum(w) ≈ DefiniteLineIntegral()*w


    ## Volterra integral equation

    d = Interval(0.0,1.0)
    V = Volterra(d)
    K = LowRankFun((x,y)->sin(y-x),d^2)
    L = I-V[K]

    @time testbandedoperator(L)

    f = Fun(exp,d)
    @test domainspace(L) == Legendre(d)
    @test rangespace(L) == Legendre(d)
    @test bandrange(V) == -1:0
    @time u = L\f
    @test norm(L*u-f) ≤ 20eps()



    ## Check DefiniteIntegral

    @time for S in (JacobiWeight(0.5,0.5,Ultraspherical(1,Segment(-2,-1))),
              JacobiWeight(0.5,0.5,Ultraspherical(1,Segment(-2,-1+2im))),
              JacobiWeight(1.5,1.5,Ultraspherical(2,Segment(-2,-1+2im))),
              JacobiWeight(-0.5,-0.5,Chebyshev(Segment(-2,-1+2im))),
              JacobiWeight(0.67,0.123,Chebyshev(Segment(-2,-1+2im))),
              JacobiWeight(0.67,0.123,Ultraspherical(1,Segment(-2,-1+2im))))
        f=Fun(S,[1.,2.,3.])
        @test DefiniteIntegral(space(f))*f ≈ sum(f)
        @test DefiniteLineIntegral(space(f))*f ≈ linesum(f)
    end


    ## Fredholm integral
    K=LowRankFun((x,y)->cos(x-y),ChebyshevInterval()^2)
    Σ=DefiniteIntegral(Chebyshev())
    testbandedbelowoperator(Σ[K])


    Σ = DefiniteIntegral(Chebyshev()); x=Fun();
    L=I+exp(x)*Σ[cos(x)]
    testbandedbelowoperator(L)






    # Piecewise space definite integral


    Γ=Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)
        z=Fun(Γ)

    S=PiecewiseSpace(map(d->isa(d,Circle) ? Fourier(d) : JacobiWeight(0.5,0.5,Ultraspherical(1,d)),components(Γ)))


    B=DefiniteLineIntegral(S)

    Random.seed!(0)
    f=Fun(S,rand(20))
    @test B*f ≈ linesum(component(f,1)) + linesum(component(f,2)) + linesum(component(f,3))



    # definite integral

    dom = ApproxFun.UnionDomain(0..1, 2..3)
    ⨍ = DefiniteLineIntegral(union(JacobiWeight.(-0.5,-0.5,ChebyshevDirichlet{1,1}.(components(dom)))...))

    x = Fun(dom)
    f = exp(x)/(sqrt(x*abs(1-x))*sqrt(abs(2-x)*abs(3-x)))
    @test sum(f) ≈ Number(⨍*f)
end
