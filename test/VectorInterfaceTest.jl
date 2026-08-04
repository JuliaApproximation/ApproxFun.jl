module VectorInterfaceTest

using ApproxFun
using Test
using VectorInterface
using VectorInterface: One, Zero

include(joinpath(@__DIR__, "testutils.jl"))

@verbose @testset "VectorInterface" begin
    f = Fun(x -> exp(x) * cos(x))
    g = Fun(x -> sin(3x))
    pts = range(-1, 1; length=10)

    @testset "scalartype" begin
        @test scalartype(f) == Float64
        @test scalartype(typeof(f)) == Float64
        @test scalartype(Fun(x -> cis(x), Fourier())) == ComplexF64
        @test scalartype(Fun(Chebyshev(), Float32[0, 1])) == Float32
    end

    @testset "zerovector" begin
        z = zerovector(f)
        @test z isa Fun
        @test space(z) == space(f)
        @test iszero(coefficients(z))
        @test scalartype(z) == Float64

        z = zerovector(f, ComplexF64)
        @test iszero(coefficients(z))
        @test scalartype(z) == ComplexF64

        h = copy(f)
        @test zerovector!(h) === h
        @test iszero(coefficients(h))

        h = copy(f)
        @test zerovector!!(h) === h
        @test iszero(coefficients(h))
    end

    @testset "scale" begin
        h = scale(f, 2.5)
        @test h(0.3) ≈ 2.5 * f(0.3)
        @test coefficients(f) == coefficients(Fun(x -> exp(x) * cos(x)))

        h = scale(f, 2im)
        @test scalartype(h) == ComplexF64
        @test h(0.3) ≈ 2im * f(0.3)

        h = copy(f)
        @test scale!(h, 3.0) === h
        @test h(0.3) ≈ 3 * f(0.3)
        @test scale!(h, One()) === h
        @test h(0.3) ≈ 3 * f(0.3)

        h = copy(f)
        @test scale!!(h, 3.0) === h
        @test h(0.3) ≈ 3 * f(0.3)

        h = copy(f)
        h2 = scale!!(h, im)
        @test h2 !== h
        @test scalartype(h2) == ComplexF64
        @test h2(0.3) ≈ im * f(0.3)

        h = copy(f)
        @test scale!(h, g, 2.0) === h
        @test ncoefficients(h) == ncoefficients(g)
        @test h(0.3) ≈ 2 * g(0.3)

        h = copy(f)
        @test scale!!(h, g, 2.0) === h
        @test h(0.3) ≈ 2 * g(0.3)

        h = copy(f)
        h2 = scale!!(h, g, im)
        @test h2 !== h
        @test h2(0.3) ≈ im * g(0.3)

        @test_throws ArgumentError scale!(copy(f), Fun(sin, 0..1), 2.0)
    end

    @testset "add" begin
        for (α, β) in ((2.0, 3.0), (2.0, One()), (One(), One()), (0.5, Zero()))
            h = add(copy(f), g, α, β)
            @test all(x -> isapprox(h(x), β * f(x) + α * g(x); atol=1e-12), pts)

            h = copy(f)
            @test add!(h, g, α, β) === h
            @test all(x -> isapprox(h(x), β * f(x) + α * g(x); atol=1e-12), pts)

            h = copy(g)
            @test add!(h, f, α, β) === h
            @test all(x -> isapprox(h(x), β * g(x) + α * f(x); atol=1e-12), pts)

            h = copy(f)
            @test add!!(h, g, α, β) === h
            @test all(x -> isapprox(h(x), β * f(x) + α * g(x); atol=1e-12), pts)
        end

        h = add(copy(f), g)
        @test all(x -> isapprox(h(x), f(x) + g(x); atol=1e-12), pts)
        h = add(copy(f), g, 2.0)
        @test all(x -> isapprox(h(x), f(x) + 2g(x); atol=1e-12), pts)
        h = add!(copy(f), g)
        @test all(x -> isapprox(h(x), f(x) + g(x); atol=1e-12), pts)

        h = copy(f)
        h2 = add!!(h, g, im, 1.0)
        @test h2 !== h
        @test scalartype(h2) == ComplexF64
        @test all(x -> isapprox(h2(x), f(x) + im * g(x); atol=1e-12), pts)

        f2 = Fun(x -> x^2, 0..1)
        g2 = Fun(x -> cospi(2x), Fourier(0..1))
        @test_throws ArgumentError add!(copy(f2), g2, 1.0, 1.0)
        h = add!!(copy(f2), g2, 2.0, 3.0)
        @test all(x -> isapprox(h(x), 3 * f2(x) + 2 * g2(x); atol=1e-10), range(0.1, 0.9; length=5))
        h = add(copy(f2), g2, 2.0, 3.0)
        @test all(x -> isapprox(h(x), 3 * f2(x) + 2 * g2(x); atol=1e-10), range(0.1, 0.9; length=5))
    end

    @testset "inner and norm" begin
        @test inner(f, g) ≈ sum(f * g)
        @test norm(f) ≈ sqrt(inner(f, f))

        u = Fun(x -> cis(x), Fourier())
        v = Fun(x -> sin(x) * cis(2x), Fourier())
        @test inner(u, v) ≈ sum(conj(u) * v)
        @test inner(v, u) ≈ conj(inner(u, v))
        @test norm(u) ≈ sqrt(real(inner(u, u)))
    end

    @testset "Gram-Schmidt orthogonalization" begin
        basis = [Fun(x -> x^k) for k in 0:3]
        ortho = typeof(f)[]
        for b in basis
            v = copy(b)
            for q in ortho
                v = add!!(v, q, -inner(q, v) / inner(q, q), One())
            end
            push!(ortho, scale!!(v, 1 / norm(v)))
        end
        for (i, p) in enumerate(ortho), (j, q) in enumerate(ortho)
            @test inner(p, q) ≈ (i == j ? 1.0 : 0.0) atol=1e-10
        end
    end
end

end # module
