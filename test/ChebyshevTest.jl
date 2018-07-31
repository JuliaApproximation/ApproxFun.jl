using ApproxFun, Test
    import ApproxFun: testspace

@testset "Chebyshev" begin
    @testset "ChebyshevInterval" begin
        @test Fun(x->4).coefficients == [4.0]
        @test Fun(4).coefficients == [4.0]

        f = Fun(ChebyshevInterval(), [1])
        @test f(0.1) == 1

        ef = Fun(exp)
        @test ef(0.1) ≈ exp(0.1)

        for d in (ChebyshevInterval(),Interval(1.,2.),Segment(1.0+im,2.0+2im))
            testspace(Chebyshev(d))
        end

        ef = Fun(exp,ChebyshevInterval())

        @test ef == -(-ef)
        @test ef == (ef-1) + 1

        ef = Fun(exp)


        @test ef == -(-ef)
        @test ef == (ef-1) + 1

        @test ef / 3 == (3 \ ef)


        cf = Fun(cos)

        ecf = Fun(x->cos(x).*exp(x))
        eocf = Fun(x->cos(x)./exp(x))

        @test ef(.5) ≈ exp(.5)
        @test ecf(.123456) ≈ cos(.123456).*exp(.123456)

        r=2 .* rand(100) .- 1

        @test maximum(abs,ef.(r)-exp.(r))<200eps()
        @test maximum(abs,ecf.(r).-cos.(r).*exp.(r))<200eps()

        @test (cf .* ef)(0.1) ≈ ecf(0.1) ≈ cos(0.1)*exp(0.1)
        @test domain(cf.*ef) ≈ domain(ecf)
        @test domain(cf.*ef) == domain(ecf)

        @test norm((ecf-cf.*ef).coefficients)<200eps()
        @test maximum(abs,(eocf-cf./ef).coefficients)<1000eps()
        @test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()

        ## Diff and cumsum
        @test norm((ef - ef').coefficients) < 10E-11
        @test norm((ef - cumsum(ef)').coefficients) < 20eps()
        @test norm((cf - cumsum(cf)').coefficients) < 20eps()
        @test sum(ef)  ≈ 2.3504023872876028
        @test norm(ef)  ≈ 1.90443178083307
    end

    @testset "Other interval" begin
        ef = Fun(exp,1..2)
        cf = Fun(cos,1..2)

        ecf = Fun(x->cos(x).*exp(x),1..2)
        eocf = Fun(x->cos(x)./exp(x),1..2)

        r=rand(100) .+ 1
        x=1.5

        @test ef(x) ≈ exp(x)

        @test maximum(abs,ef.(r)-exp.(r))<400eps()
        @test maximum(abs,ecf.(r).-cos.(r).*exp.(r))<100eps()
        @test norm((ecf-cf.*ef).coefficients)<500eps()
        @test maximum(abs,(eocf-cf./ef).coefficients)<1000eps()
        @test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()

        ## Diff and cumsum
        @test norm((ef - ef').coefficients)<10E-11
        @test norm((ef - cumsum(ef)').coefficients) < 10eps()
        @test norm((cf - cumsum(cf)').coefficients) < 10eps()
        @test sum(ef) ≈ 4.670774270471604
        @test norm(ef) ≈ 4.858451087240335
    end

    @testset "Roots" begin
        f=Fun(x->sin(10(x-.1)))
        @test norm(f.(roots(f)))< 1000eps()

        @test_throws ArgumentError roots(Fun(zero))
        @test_throws ArgumentError roots(Fun(Chebyshev(),Float64[]))
    end

    @testset "Aliasing" begin
        f=Fun(x->cos(50acos(x)))
        @test norm(f.coefficients-Matrix(I,ncoefficients(f),ncoefficients(f))[:,51])<100eps()
    end

    @testset "Integer values" begin
        @test Fun(x->2,10)(.1) ≈ 2
        @test Fun(x->2)(.1) ≈ 2

        @test Fun(Chebyshev,Float64[]).([0.,1.]) ≈ [0.,0.]
        @test Fun(Chebyshev,[])(0.) ≈ 0.
    end

    @testset "large intervals #121" begin
        x = Fun(identity,0..10)
        f = sin(x^2)
        g = cos(x)
        @test f(.1) ≈ sin(.1^2)

        x = Fun(identity,0..100)
        f = sin(x^2)
        @test ≈(f(.1),sin(.1^2);atol=1E-12)
    end

    @testset "Reverse" begin
        f=Fun(exp)
        @test ApproxFun.default_Fun(f, Chebyshev(Segment(1 , -1)), ncoefficients(f))(0.1) ≈ exp(0.1)
        @test Fun(f,Chebyshev(Segment(1,-1)))(0.1) ≈ f(0.1)
    end

    @testset "minimum/maximum" begin
        x=Fun()
        @test minimum(x) == -1
        @test maximum(x) == 1
    end

    @testset "Do not overresolve #7" begin
        @test ncoefficients(Fun(x->sin(400*pi*x),-1..1)) ≤ 1400
    end

    @testset "Bug from Trogdon" begin
        δ = .03 # should be less than 0.03
        @test 0.0 ∈ Domain(1-8*sqrt(δ)..1+8*sqrt(δ))
        @test 0.00001 ∈ Domain(1-8*sqrt(δ)..1+8*sqrt(δ))

        ϕfun = Fun(x -> 1/sqrt(2*pi*δ)*exp(-abs2.(x-1)/(2*δ)), 1-8sqrt(δ)..1+8sqrt(δ))
        ϕfun(0.00001) ≈ 1/sqrt(2*pi*δ)*exp(-abs2.(0.00001-1)/(2*δ))

        iϕfun = 1-cumsum(ϕfun)
        @test iϕfun(0.00001) ≈ 1
    end

    @testset "Large scaling" begin
        w = Fun(x -> 1e5/(x*x+1), 283.72074879785936 .. 335.0101119042838)
        @test w(leftendpoint(domain(w))) ≈ 1e5/(leftendpoint(domain(w))^2+1)
    end
end
