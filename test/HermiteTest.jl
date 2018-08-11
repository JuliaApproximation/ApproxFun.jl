using ApproxFun, Compat.Test
    import ApproxFun: testbandedoperator

@testset "Hermite and GaussWeight" begin
    @testset "Evaluation" begin
        f = Fun(x-> x + x^2, Hermite())
        @test f(1.0) ≈ 2.0
        @test values(f) ≈ points(f) + points(f).^2

        w = Fun(GaussWeight(), [1.0])
        @test w(0.1) ≈ exp(-0.1^2)
        w = Fun(GaussWeight(), Float64[])
        @test w(1) == 0

        w = Fun(GaussWeight(), Float64[1.0])
        f = Fun(x-> 1 + x + x^2, Hermite()) * w
        @test f(0.1) ≈ exp(-0.1^2) * (1+0.1+0.1^2)

        w = Fun(GaussWeight(Hermite(2), 0), [1.0,2.0,3.0]);
        w̃ = Fun(Hermite(2), [1.0,2.0,3.0]);
        @test w(0.1) == w̃(0.1)

        w = Fun(GaussWeight(Hermite(2), 2), Float64[1.0])
        f = Fun(x-> 1 + x + x^2, Hermite()) * w
        @test f(0.1) ≈ exp(-0.1^2) * (1+0.1+0.1^2)

        L = 1.3; x = 1.2;
        H₀ = Fun(Hermite(), [1.0])
        H̃₀ = Fun(Hermite(L), [1.0])
        @test H̃₀(x) ≈ H₀(sqrt(L) * x)

        L = 1.3; x = 1.2;
        H₁ = Fun(Hermite(), [0.0,1.0])
        H̃₁ = Fun(Hermite(L), [0.0,1.0])
        @test H̃₁(x) ≈ H₁(sqrt(L) * x)
    end


    @testset "Derivative" begin
        D = Derivative(Hermite())
        testbandedoperator(D)

        f = Fun( x-> x + x^2, Hermite())
        g = D * f
        @test g(1.) ≈ 3.
    end


    @testset "Integration" begin
        @test_throws ArgumentError integrate(Fun(GaussWeight(Hermite(2),1), [0.0,1.0]))

        w = Fun(GaussWeight(Hermite(2), 0), [1.0,2.0,3.0])
        g = integrate(w)
        g̃ = Fun(Hermite(2), [0.0, 0.5, 0.5, 0.5])
        @test g(0.1) == g̃(0.1)

        w = Fun(GaussWeight(), Float64[])
        g = integrate(w)
        @test g(0.1) == 0.0

        w = Fun(GaussWeight(), [1.0])
        g = integrate(w)
        @test_skip w̃ = Fun(w, -7..7)
        w̃ = Fun( x-> w(x), -7..7)
        g̃ = cumsum(w̃)
        @test g(3) - g(-7) ≈ g̃(3)

        w = Fun(GaussWeight(), Float64[1.0])
        g = integrate(w)
        @test_skip w̃ = Fun(w, -7..7)
        w̃ = Fun(x -> w(x), -7..7)
        g̃ = cumsum(w̃)
        @test g(3) - g(-7) ≈ g̃(3)

        w = Fun(GaussWeight(Hermite(2), 2), Float64[1.0])
        g = integrate(w)
        @test_skip w̃ = Fun(w, -7..7)
        w̃ = Fun(x -> w(x), -7..7)
        g̃ = cumsum(w̃)
        @test g(3) - g(-7) ≈ g̃(3)

        w = Fun(GaussWeight(), Float64[0.0, 1.0])
        g = integrate(w)
        @test_skip w̃ = Fun(w, -7..7)
        w̃ = Fun(x -> w(x), -7..7)
        g̃ = cumsum(w̃)
        @test g(3) - g(-7) ≈ g̃(3)

        w = Fun(GaussWeight(Hermite(2), 2), Float64[0.0, 1.0])
        g = integrate(w)
        @test_skip w̃ = Fun(w, -7..7)
        w̃ = Fun(x -> w(x), -7..7)
        g̃ = cumsum(w̃)
        @test g(3) - g(-7) ≈ g̃(3)
    end
end
