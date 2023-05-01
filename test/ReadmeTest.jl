module ReadmeTests

using ApproxFun
using SpecialFunctions
using LinearAlgebra
using Test

const EXAMPLES_DIR = joinpath(dirname(dirname(pathof(ApproxFun))), "examples")

include(joinpath(@__DIR__, "testutils.jl"))

function get_included_files(filename)
    v = [l for l in eachline(joinpath(EXAMPLES_DIR, filename)) if contains(l, "include")]
    strip.(getindex.(split.(getindex.(split.(v, "include("), 2), ")"), 1), '\"')
end

@verbose @testset "Readme" begin
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

    @testset "NonlinearBVP" begin
        for f in get_included_files("NonlinearBVP.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end

    @testset "ODE" begin
        for f in get_included_files("ODE.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end
    @testset "System of Equations" begin
        for f in get_included_files("System.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end

    @testset "Periodic" begin
        for f in get_included_files("Periodic.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end
    @testset "Sampling" begin
        for f in get_included_files("Sampling.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end
    @testset "Eigenvalue" begin
        for f in get_included_files("Eigenvalue.jl")
            include(joinpath(EXAMPLES_DIR, f))
        end
    end

    # this is broken on v1.6 (on BandedArrays < v0.16.16), so we only test on higher versions
    if VERSION >= v"1.7"
        @testset "PDE" begin
            for f in get_included_files("PDE.jl")
                include(joinpath(EXAMPLES_DIR, f))
            end
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

end # module
