using ApproxFun, SpecialFunctions, Test



@testset "broadcast" begin
    @testset "Special functions" begin
        x = Fun()
        @test exp(x) ≈ exp.(x) atol=10eps()
        f = Fun(exp)
        @test exp.(f) ≈ exp(f)
        @test besselj.(1,f) ≈ besselj(1,f)

        @test atan.(f,1)(0.1) ≈ atan(f(0.1),1)
        @test atan.(f,f)(0.1) ≈ atan(f(0.1),f(0.1))
    end
    
    @testset  "triplet with mixed scalar-fun" begin
        f = Fun(exp)
        @test ((a,b,c) -> exp(a +b+c)).(f,2.0,f) ≈ exp(2f + 2.0)
    end

    @testset "Arrays" begin
        x = Fun()
        @test exp.( x .+ [1,2,3]) isa Vector
        @test exp.(x .+ [1,2]) ≈ [exp(x+1),exp(x+2)]

        @test maximum.( [x;x] .+ [x;x]) isa Fun
        @test maximum.( [x;x] .+ [x;x]) ≈ 2x

        @test exp.([x,x] .+ 2 .+ [x;x]) isa Vector
        @test exp.( [x,x] .+ 2 .+ [x;x]) ≈ [exp(2x+2),exp(2x+2)]


        # since array broadcasting takes presidence, the following
        # does not try to call a constructor:
        @test ApproxFun.tocanonical.(x,[0.1,0.2]) ≈ [0.1,0.2]
    end

    @testset "broadcast!" begin
        x = Fun()
        f = Fun()
        f .= exp.(x)
        @test ≈(exp(x),f;atol=10eps())

        f = Fun(Ultraspherical(1))
        f .= exp.(x)
        @test f(0.1) ≈ exp(0.1)

        f = Fun()
        @test_throws ArgumentError (f .= Fun(Line()))
    end

    @testset "norm" begin
        f = Fun(x ->[exp(x),sin(x)])
        F = [Fun(exp),Fun(sin)]
        @test norm.(f) ≈ Fun(x ->norm([exp(x),sin(x)]))
        @test norm.(F) ≈ [norm(Fun(exp)),norm(Fun(sin))] # isa Vector{Float64}

        @test norm.(f,[1,2]) ≈ norm.(F,[1,2]) ≈ [norm(f[1],1),norm(f[2],2)] # isa Vector{Float64}
    end
end
