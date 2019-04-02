using ApproxFun, Test

@time include("ReadmeTest.jl")
@time include("ExtrasTest.jl")
@time include("NumberTypeTest.jl")
@time include("FractionalTest.jl")

@testset "Chebyshev and Fourier" begin
    @test norm(Fun(x->Fun(cos,Fourier(-π .. π),20)(x),20)-Fun(cos,20)) <100eps()
    @test norm(Fun(x->Fun(cos,Fourier(-π .. π))(x))-Fun(cos)) <100eps()
    @test norm(Fun(x->Fun(cos,Laurent)(x))-Fun(cos)) <100eps()
end


@testset "Piecewise space definite integral" begin
    Γ=Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)
        z=Fun(Γ)

    S=PiecewiseSpace(map(d->isa(d,Circle) ? Fourier(d) : JacobiWeight(0.5,0.5,Ultraspherical(1,d)),components(Γ)))


    B=DefiniteLineIntegral(S)

    Random.seed!(0)
    f=Fun(S,rand(20))
    @test B*f ≈ linesum(component(f,1)) + linesum(component(f,2)) + linesum(component(f,3))
end

@testset "Extending function" begin
    Γ = Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)

    @test isempty(component(Γ,1)\component(Γ,1))
    @test Γ \ component(Γ,1) == component(Γ,2) ∪ component(Γ,3)

    @test norm(Fun(ones(component(Γ,1)),Γ) - Fun(x->x ∈ component(Γ,1) ? 1.0 : 0.0,Γ)) == 0
end

@testset "periodic x interval" begin
    dθ=PeriodicSegment(-2.,2.);dt=Interval(0,1.)
    d=dθ×dt
    Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
    u0=Fun(θ->exp(-20θ^2),dθ)

    ε = 0.1
    @time u=\([I⊗ldirichlet(dt);Dt-ε*Dθ^2-Dθ],[u0;0.];tolerance=1E-4)
    @test ≈(u(0.1,0.2),0.3103472600253807;atol=1E-2)

    A=Dt+Dθ
    testbandedblockbandedoperator(A)

    @time u=\([I⊗ldirichlet(dt);Dt+Dθ],[u0;0.0];tolerance=1E-6)
    @test ≈(u(0.2,0.1),u0(0.1);atol=1E-6)
end


@testset "Laplace in a strip" begin
    d=PeriodicSegment() × ChebyshevInterval()
    g=Fun((x,y)->real(cos(x+im*y)),∂(d))  # boundary data


    @test g(0.1,1.0) ≈ real(cos(0.1+im))
    @test g(0.1,-1.0) ≈ real(cos(0.1-im))
    v=[g;0]
    @test v(0.1,-1) ≈ [real(cos(0.1-im));0]

    A=[Dirichlet(d);Laplacian(d)]
    a = space(v)
    b = rangespace(A)

    @test Fun(component(v[1],1), component(b[1],1))(0.1,-1.0) ≈ v(0.1,-1.0)[1]
    @test Fun(component(v[1],2), component(b[1],2))(0.1,-1.0) ≈ v(0.1,-1.0)[1]
    @test ApproxFun.default_Fun(v[1] , b[1])(0.1,1.0) ≈ v(0.1,1.0)[1]
end

@testset "Bug in Multiplication" begin
    dom = Interval(0.001, 1) × PeriodicSegment(-pi, pi)

    @test blocklengths(Space(dom)) == 2:2:∞

    r,r2 = Fun((r,t) -> [r;r^2], dom)

    @test r(0.1,0.2) ≈ 0.1
    @test r2(0.1,0.2) ≈ 0.1^2

    sp = Space(dom)
    Dr = Derivative(sp, [1,0])
    @test ApproxFun.blockbandwidths(Dr) == (-1,1)
    @test ApproxFun.subblockbandwidths(Dr)  == (1,3)

    Dθ = Derivative(sp, [0,1])
    Mr = Multiplication(Fun( (r, θ) -> r, sp ), sp)
    rDr = Mr * Dr

    testbandedblockbandedoperator(rDr)
end

@testset "Periodic x Interval" begin
    d=PeriodicSegment() × ChebyshevInterval()

    u_ex=Fun((x,y)->real(cos(x+im*y)),d)
    @test u_ex(1.0,0.1) ≈ real(cos(1.0+im*0.1)) atol=10eps()

    B=Dirichlet(Space(d))

    @test B.order == 0  # tests stupid bug
    g=Fun((x,y)->real(cos(x+im*y)),rangespace(B))  # boundary data

    @test norm((B*u_ex-g).coefficients) < 100eps()

    testbandedblockbandedoperator(Laplacian(d))

    @time u=[B;Laplacian(d)]\[g;0.]

    @test u(.1,.2) ≈ real(cos(.1+.2im))
end

@testset "Mix Fourier-Chebyshev (#602)" begin
    s = Chebyshev(-π..π)
    a = Fun(t-> 1+sin(cos(2t)), s)
    L = Derivative() + a
    f = Fun(t->exp(sin(10t)), s)
    B = periodic(s,0)
    @time uChebyshev = [B;L] \ [0.;f]

    s = Fourier(-π..π)
    a = Fun(t-> 1+sin(cos(2t)), s)
    L = Derivative() + a
    f = Fun(t->exp(sin(10t)), s)
    @time uFourier = L\f

    @test norm(uFourier-uChebyshev) ≤ 100eps()
end

@testset "Conversion" begin
    f=Fun(t->[cos(t) 0;sin(t) 1],-π..π)
    g=Fun(f,Space(PeriodicSegment(-π,π)))
    @test g(.1) ≈ f(.1)
end

@testset "definite integral" begin
    Σ = DefiniteIntegral()

    f1 = Fun(t->cos(cos(t)),-π..π)
    f = Fun(t->cos(cos(t)),Laurent(-π..π))

    @test sum(f1) ≈ Σ*f
end


@testset "Sampling" begin
    ff=(x,y)->(x-y)^2*exp(-x^2/2-y^2/2)

    f=Fun(ff,Domain(-4..4)^2)
    r=ApproxFun.sample(f,5000)


    #We can compare the histogram to the 1-point correlation
    g=sum(f,1)/sum(f)
    @test  g(0.1) ≈ 0.2004758624973169

    # check bisection inv
    f = Fun(x -> exp(-x^2/2),-5..5)
    g = cumsum(f)
    @test g(ApproxFun.bisectioninv(g,0.5)) ≈ 0.5
end

@testset "piecewise sample (#635)" begin
    f = abs(Fun(sin, -5..5))
    @test integrate(f)(-4.0) ≈ -(cos(-4.0) - cos(-5.0))
    @test -(cos(-π) - cos(-5.0)) + cos(-3.0) - cos(-π) ≈ integrate(f)(-3.0)
    r = ApproxFun.sample(f,10)
    @test maximum(r) ≤ 5
    @test minimum(r) ≥ -5
end

f=Fun(exp)
x=sample(f,100000)
x=sample(f,100000)
@time x=sample(f,100000)
println("Sample: Time should be ~0.13")
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds

