module FractionalTest

using ApproxFun
using Test
using LinearAlgebra
using SpecialFunctions
using ApproxFunBaseTest: testfunctional, testbandedoperator
include(joinpath(@__DIR__, "testutils.jl"))

@verbose @testset "Fractional" begin
    @testset "Jupyter example" begin
        S = Legendre() ⊕ JacobiWeight(0.5,0.,Ultraspherical(1))
        @time Q½ = LeftIntegral(S,0.5)
        @time testbandedoperator(Q½)
        @time I+Q½
        @time testbandedoperator(I+Q½)
        @time y = (I+Q½)\1
        @test values(y)[1] ≈ 0.33627096683893143

        S = Legendre()⊕JacobiWeight(0.5,0.,Ultraspherical(1))
        testfunctional(rdirichlet(S))
        B=rdirichlet(S)
        D=LeftDerivative(1.5) : S
        testbandedoperator(D)
    end

    ## Avazzadev et al

    @testset "Example 1" begin
        x=Fun(0..1)
        Q=gamma(0.5)*LeftIntegral(0.5)
        @time f=(2/105*sqrt(x)*(105-56x^2+48x^3))
        u=Q\f
        @test norm(u-(x^3-x^2+1))<100eps()
    end


    @testset "Example 2" begin
        x=Fun(0..1)
        Q=gamma(0.5)*LeftIntegral(0.5)
        u=Q\(exp(x)-1)
        @time @test norm(u-exp(x)*erf(sqrt(x))/sqrt(π)) < 100eps() # 5.0036177384681187e-14
    end


    @testset "Example 3" begin
        x=Fun(0..1)
        Q=gamma(1/5)*LeftIntegral(1/5)
        u=Q\(x+1)
        @test norm(u-(1+1.25x)*sin(0.8π)/(π*x^(1/5))) < 10eps()
    end

    @testset "Example 4" begin
        x=Fun(0..1)
        Q=gamma(1-1/3)*LeftIntegral(1-1/3)
        u=Q\x^(7/6)
        @test norm(u-7*gamma(1/6)/(18*sqrt(π)*gamma(2/3))*sqrt(x)) < 100eps()
    end

    @testset "Example 5" begin
        d=Interval(0,1)
        x=Fun(d)
        f=x+4/3*x^(3/2)
        S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
        Q=gamma(.5)*LeftIntegral(S,.5)
        @time @test sum(f/sqrt(1-x)) ≈ last(Q*f)

        L=I+Q
        @test last(L.ops[2]*f) ≈ last(Q*f)
        @test last(L*f) ≈ last(f)+last(Q*f)

        @time u=L\f
        @test norm(u-x)  < 10eps()
    end

    @testset "Example 6" begin
        d=Interval(0,1)
        x=Fun(d)
        @time f=x^2+16/15*x^(5/2)
        S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
        Q=gamma(.5)*LeftIntegral(S,.5)
        L=I+Q
        @time u=L\f
        @test norm(u-x^2) < 10eps()
    end

    @testset "Example 7" begin
        d=Interval(0.,1.)
        x=Fun(d)
        f=2sqrt(x)
        S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
        Q=gamma(.5)*LeftIntegral(S,.5)
        @time L=I+Q
        u=L\f

        @time @test norm(1-exp(π*x)*erfc(sqrt(π*x))-u) < 100eps()
    end


    @testset "Example 8" begin
        d=Interval(0.,1.)
        x=Fun(d)
        @time f=1/(x+1)+2*Fun(x->asinh(sqrt(x))/sqrt(1+x),JacobiWeight(.5,0.,d))
        S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
        Q=gamma(.5)*LeftIntegral(S,.5)
        L=I+Q
        u=L\f
        @test norm((u-1/(x+1)).coefficients) < 1000eps()   # 1.2011889731154679e-14
    end



    @testset "Test for bug" begin
        QL = LeftIntegral(0.5)	: Legendre() →	JacobiWeight(0.5,0.,Ultraspherical(1))
        QU = LeftIntegral(0.5)	: JacobiWeight(0.5,0.,Ultraspherical(1)) → Legendre()

        λ=0.25
        @time L=[λ*I QU; QL λ*I]
        @test L[2,5] ≈ 0.
    end

    @testset "paper examples" begin
        S=Legendre()⊕JacobiWeight(0.5,0.,Ultraspherical(1))
        Q½=LeftIntegral(S,0.5)

        y=(I+Q½)\1

        x=Fun()
        @time @test norm(exp(1+x)*erfc(sqrt(1+x))-y) < 100eps()

        S=Legendre()⊕JacobiWeight(0.5,0.,Ultraspherical(1))
        x=Fun()
        Q½=LeftIntegral(S,0.5)

        @time y=(I+exp(-(1+x)/2)*Q½[exp((1+x)/2)])\exp(-(1+x)/2)

        @test norm(y-exp((1+x)/2)*erfc(sqrt(1+x))) < 100eps()
    end

    @testset "type inference" begin
        if VERSION >= v"1.8"
            @inferred (() -> LeftIntegral(Jacobi(-1,0)))()
            @inferred (() -> LeftIntegral(Legendre()))()
            # @inferred (() -> LeftIntegral(Chebyshev()))()
            @inferred (() -> LeftIntegral(Ultraspherical(0.5)))()
            @inferred (() -> LeftIntegral(JacobiWeight(0,0,Legendre())))()
            @inferred (() -> LeftIntegral(JacobiWeight(-0.5,0,Chebyshev())))()
            @inferred (() -> RightIntegral(Jacobi(-1,0)))()
            @inferred (() -> RightIntegral(Legendre()))()
            @inferred (() -> RightIntegral(JacobiWeight(0,0,Legendre())))()
            @inferred (() -> RightIntegral(JacobiWeight(0,-0.5,Chebyshev())))()
            @inferred (() -> LeftIntegral(Legendre() ⊕ JacobiWeight(0.5,0.,Ultraspherical(1)),0.5))()
        end
    end
end

end # module
