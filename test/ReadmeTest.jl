using ApproxFun, SpecialFunctions, LinearAlgebra, Test


@testset "Readme" begin
    @testset "Calculus and algebra" begin
        x = Fun(identity,0..10)
        f = sin(x^2)
        g = cos(x)

        @test ≈(f(.1),sin(.1^2);atol=1000eps())

        h = f + g^2
        r = roots(h)
        rp = roots(differentiate(h))

        @test norm(h.(r))<1000eps()
        @test norm(h'.(rp))<100000eps()
    end



    @testset "Differentiation and Integration" begin
        f = Fun(x->exp(x),-1..1)
        @test norm(f-f')<1000eps()

        g = cumsum(f)
        g = g + f(-1)
        @test norm(f-g)<100eps()



        x = Fun(identity)
        f = exp(x)
        g = f/sqrt(1-x^2)


        space(f),domain(f)
        space(g),domain(g)


        f = Fun(x->cospi(5x))
        g = abs(f)
        space(f)
        space(g)

        x = Fun()
        f = erf(x)
        g = besselj(3,exp(f))
        h = airyai(10asin(f)+2g)
    end

    @testset "ODE" begin
        a,b = -1000,200
        x = Fun(identity, a..b)
        d = domain(x)
        D = Derivative(d)
        B = Dirichlet(d)
        L = D^2 - x
        u = [B;L] \ [[airyai(a),airyai(b)],0]

        @test ≈(u(0.),airyai(0.);atol=10000eps())


        ##) Nonlinear BVPs
        x=Fun()
        u0=0.0x

        N=u->[u(-1.)-1.,u(1.)+0.5,0.001u''+6*(1-x^2)*u'+u^2-1.]
        u=newton(N,u0)

        @test norm(N(u)[end]) ≤ 1000eps()
    end


    @testset "Periodic" begin
        f = Fun(cos,Fourier(-π..π))
        @test norm(differentiate(f) + Fun(sin,Fourier(-π..π))) < 100eps()

        s = Chebyshev(-π..π)
        a = Fun(t-> 1+sin(cos(2t)),s)
        L = Derivative() + a
        f = Fun(t->exp(sin(10t)),s)
        B = periodic(s,0)
        uChebyshev = [B;L]\[0.,f]

        s = Fourier(-π..π)
        a = Fun(t-> 1+sin(cos(2t)),s)
        L = Derivative() + a
        f = Fun(t->exp(sin(10t)),s)
        uFourier = L\f


        @test uChebyshev(0.) ≈ uFourier(0.)
    end

    @testset "Sampling" begin
        ## Sampling
        f = abs(Fun(sin,-5..5))
        x = ApproxFun.sample(f,10)
    end


    @testset "PDE" begin
        d = ChebyshevInterval()^2                            # Defines a rectangle

        # @time u = \([Dirichlet(d);Laplacian(d)+100I],
        #                     [ones(∂(d));0.];tolerance=1E-10)      # First four entries of rhs are
        #

        # this is broken on v1.6 (on BandedArrays < v0.16.16), so we skip on errors
        if VERSION >= v"1.7"
            QR = qr([Dirichlet(d);Laplacian()+100I])
            @time ApproxFun.resizedata!(QR,:,4000)
            @time u = \(QR, [ones(∂(d)); 0.]; tolerance=1E-7)
            @test u(0.1,1.) ≈ 1.0
            @test ≈(u(0.1,0.2),-0.02768276827514463;atol=1E-8)
        end
    end

    # @testset "BigFloat" begin
    @test_skip begin
        setprecision(1000) do
            d=BigFloat(0)..BigFloat(1)
            D=Derivative(d)
            u=[ldirichlet();D-I]\[1;0]
            @test u(1) ≈ exp(BigFloat(1))
        end
    end
end
