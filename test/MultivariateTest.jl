using ApproxFun, LinearAlgebra, SpecialFunctions, Test
    import ApproxFun: testbandedblockbandedoperator, testraggedbelowoperator, factor, Block, cfstype,
                        blocklengths, block, tensorizer, Vec, ArraySpace, ∞,
                        testblockbandedoperator


@testset "Multivariate" begin
    @testset "Square" begin
        S = Space(ChebyshevInterval()^2)
        @test @inferred(blocklengths(S)) ≡ Base.OneTo(∞)

        @test block(tensorizer(S), 1) == Block(1)

        @time for k=0:5,j=0:5
            ff=(x,y)->cos(k*acos(x))*cos(j*acos(y))
            f=Fun(ff,ChebyshevInterval()^2)
            @test f(0.1,0.2) ≈ ff(0.1,0.2)
        end

        @time for k=0:5,j=0:5
            ff = (x,y)->cos(k*acos(x/2))*cos(j*acos(y/2))
            f=Fun(ff,Interval(-2,2)^2)
            @test f(0.1,0.2) ≈ ff(0.1,0.2)
        end


        @time for k=0:5,j=0:5
            ff=(x,y)->cos(k*acos(x-1))*cos(j*acos(y-1))
            f=Fun(ff,Interval(0,2)^2)
            @test f(0.1,0.2) ≈ ff(0.1,0.2)
        end

        ## Try constructor variants

        ff=(x,y)->exp(-10(x+.2)^2-20(y-.1)^2)*cos(x*y)
        gg=x->exp(-10(x[1]+.2)^2-20(x[1]-.1)^2)*cos(x[1]*x[2])
        f=Fun(ff,ChebyshevInterval()^2,10000)
        @test f(0.,0.) ≈ ff(0.,0.)


        f=Fun(gg,ChebyshevInterval()^2,10000)
        @test f(0.,0.) ≈ ff(0.,0.)

        f=Fun(ff,ChebyshevInterval()^2)
        @test f(0.,0.) ≈ ff(0.,0.)
        f=Fun(gg,ChebyshevInterval()^2)
        @test f(0.,0.) ≈ ff(0.,0.)

        f=Fun(ff)
        @test f(0.,0.) ≈ ff(0.,0.)
        f=Fun(gg)
        @test f(0.,0.) ≈ ff(0.,0.)


        # Fun +-* constant
        f=Fun((x,y)->exp(x)*cos(y))

        @test f(0.1,0.2)+2 ≈ (f+2)(0.1,0.2)
        @test f(0.1,0.2)-2 ≈ (f-2)(0.1,0.2)
        @test f(0.1,0.2)*2 ≈ (f*2)(0.1,0.2)
    end

    @testset "LowRankFun" begin
        @time F = LowRankFun((x,y)->besselj0(10(y-x)),Chebyshev(),Chebyshev())

        @test F(.123,.456) ≈ besselj0(10(.456-.123))

        @time G = LowRankFun((x,y)->besselj0(10(y-x));method=:Cholesky)

        @test G(.357,.246) ≈ besselj0(10(.246-.357))

        # test "fast" grid evaluation of LowRankFun
        f = LowRankFun((x,y) -> exp(x) * cos(y)); n = 1000
        x = range(-1, stop=1, length=n); y = range(-1, stop=1, length=n)
        X = x * fill(1.0,1,n); Y = fill(1.0, n) * y'
        @time v1 = f.(X, Y);
        @time v2 = f.(x, y');
        @test v1 ≈ v2
    end


    @testset  "Vec segment" begin
        d = Segment(Vec(0.,0.) , Vec(1.,1.))
        x = Fun()
        @test (ApproxFun.complexlength(d)*x/2)(0.1)  ≈ (d.b - d.a)*0.1/2
        @test ApproxFun.fromcanonical(d,x)(0.1) ≈ (d.b+d.a)/2 + (d.b - d.a)*0.1/2

        x,y = Fun(Segment(Vec(0.,0.) , Vec(2.,1.)))
        @test x(0.2,0.1) ≈ 0.2
        @test y(0.2,0.1) ≈ 0.1

        d=Segment((0.,0.),(1.,1.))
        f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d)
        @test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))
        @test f(Vec(0.5,0.5)) ≈ exp(-0.5-2cos(0.5))

        f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d,20)
        @test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))

        f=Fun((x,y)->exp(-x-2cos(y)),d)
        @test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))

        f=Fun((x,y)->exp(-x-2cos(y)),d,20)
        @test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))
    end

    @testset "Multivariate calculus" begin
        ## Sum
        ff = (x,y) -> (x-y)^2*exp(-x^2/2-y^2/2)
        f=Fun(ff, (-4..4)^2)
        @test f(0.1,0.2) ≈ ff(0.1,0.2)

        @test sum(f,1)(0.1) ≈ 2.5162377980828357
        f=LowRankFun(f)
        @test evaluate(f.A,0.1) ≈ map(f->f(0.1),f.A)
    end


    @testset "KroneckerOperator" begin
        Mx = Multiplication(Fun(cos),Chebyshev())
        My = Multiplication(Fun(sin),Chebyshev())
        K = Mx⊗My

        @test ApproxFun.BandedBlockBandedMatrix(view(K,1:10,1:10)) ≈ [K[k,j] for k=1:10,j=1:10]
        C = Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical(1)⊗Ultraspherical(1))
        @test C[1:100,1:100] ≈ Float64[C[k,j] for k=1:100,j=1:100]
    end

    @time @testset "Partial derivative operators" begin
        d = Space(0..1) * Space(0..2)
        Dx = Derivative(d, [1,0])
        testbandedblockbandedoperator(Dx)
        f = Fun((x,y) -> sin(x) * cos(y), d)
        fx = Fun((x,y) -> cos(x) * cos(y), d)
        @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
        Dy = Derivative(d, [0,1])
        testbandedblockbandedoperator(Dy)
        fy = Fun((x,y) -> -sin(x) * sin(y), d)
        @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
        L = Dx + Dy
        testbandedblockbandedoperator(L)
        @test (L*f)(0.2,0.3) ≈ (fx(0.2,0.3)+fy(0.2,0.3))

        B=ldirichlet(factor(d,1))⊗ldirichlet(factor(d,2))
        @test abs(Number(B*f)-f(0.,0.)) ≤ 10eps()

        B=Evaluation(factor(d,1),0.1)⊗ldirichlet(factor(d,2))
        @test Number(B*f) ≈ f(0.1,0.)

        B=Evaluation(factor(d,1),0.1)⊗Evaluation(factor(d,2),0.3)
        @test Number(B*f) ≈ f(0.1,0.3)

        B=Evaluation(d,(0.1,0.3))
        @test Number(B*f) ≈ f(0.1,0.3)
    end

    @time @testset "x,y constructors" begin
        d=ChebyshevInterval()^2

        sp = ArraySpace(d,2)
        @test blocklengths(sp) == 2:2:∞
        @test block(ArraySpace(d,2),1) == Block(1)

        x,y=Fun(d)
        @test x(0.1,0.2) ≈ 0.1
        @test y(0.1,0.2) ≈ 0.2

        x,y=Fun(identity, d, 20)
        @test x(0.1,0.2) ≈ 0.1
        @test y(0.1,0.2) ≈ 0.2


        # Boundary

        x,y=Fun(identity, ∂(d), 20)
        @test x(0.1,1.0) ≈ 0.1
        @test y(1.0,0.2) ≈ 0.2


        x,y=Fun(identity, ∂(d))
        @test x(0.1,1.0) ≈ 0.1
        @test y(1.0,0.2) ≈ 0.2


        x,y=Fun(∂(d))
        @test x(0.1,1.0) ≈ 0.1
        @test y(1.0,0.2) ≈ 0.2
    end

    @testset "conversion between" begin
        dx = dy = ChebyshevInterval()
        d = dx × dy
        x,y=Fun(d)
        @test x(0.1,0.2) ≈ 0.1
        @test y(0.1,0.2) ≈ 0.2

        x,y = Fun(∂(d))
        x,y = components(x),components(y)

        g = [real(exp(x[1]-1im));0.0y[2];real(exp(x[3]+1im));real(exp(-1+1im*y[4]))]
        B = [ Operator(I,dx)⊗ldirichlet(dy);
             ldirichlet(dx)⊗Operator(I,dy);
             Operator(I,dx)⊗rdirichlet(dy);
             rneumann(dx)⊗Operator(I,dy)    ]


        @test Fun(g[1],rangespace(B)[1])(-0.1,-1.0) ≈ g[1](-0.1,-1.0)
        @test Fun(g[3],rangespace(B)[3])(-0.1,1.0)  ≈ g[3](-0.1,1.0)


        A = [B; Laplacian()]

        @test cfstype([g;0.0]) == Float64
        g2 = Fun([g;0.0],rangespace(A))
        @test cfstype(g2) == Float64

        @test g2[1](-0.1,-1.0) ≈ g[1](-0.1,-1.0)
        @test g2[3](-0.1,1.0)  ≈ g[3](-0.1,1.0)

        S=WeightedJacobi(1,1)^2
        L=Laplacian(S)
        testbandedblockbandedoperator(L)
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

    @testset "Cheby * Interval" begin
        d = ChebyshevInterval()^2
        x,y = Fun(∂(d))

        @test ApproxFun.rangetype(Space(∂(d))) == Float64
        @test ApproxFun.rangetype(space(y)) == Float64

        @test (im*y)(1.0,0.1) ≈ 0.1im
        @test (x+im*y)(1.0,0.1) ≈ 1+0.1im

        @test exp(x+im*y)(1.0,0.1) ≈ exp(1.0+0.1im)
    end


    @testset "Taylor()^2, checks bug in type of plan_transform" begin
        f = Fun((x,y)->exp((x-0.1)*cos(y-0.2)),Taylor()^2)
        @test f(0.2,0.3) ≈ exp(0.1*cos(0.1))
    end

    @testset "DefiniteIntegral" begin
        f = Fun((x,y) -> exp(-x*cos(y)))
        @test Number(DefiniteIntegral()*f) ≈ sum(f)
    end

    @testset "Piecewise Tensor" begin
        a = Fun(0..1) + Fun(2..3)
        f = a ⊗ a
        @test f(0.1,0.2) ≈ 0.1*0.2
        @test f(1.1,0.2) ≈ 0
        @test f(2.1,0.2) ≈ 2.1*0.2

        @test component(space(f),1,1) == Chebyshev(0..1)^2
        @test component(space(f),1,2) == Chebyshev(0..1)*Chebyshev(2..3)
        @test component(space(f),2,1) == Chebyshev(2..3)*Chebyshev(0..1)
        @test component(space(f),2,2) == Chebyshev(2..3)^2
    end

    @testset "Bug in chop of ProductFun" begin
        u = Fun(Chebyshev()^2,[0.0,0.0])
        @test coefficients(chop(ProductFun(u),10eps())) == zeros(0,1)


        d= (-1..1)^2
        B=[Dirichlet(factor(d,1))⊗I;I⊗ldirichlet(factor(d,2));I⊗rneumann(factor(d,2))]
        Δ=Laplacian(d)

        rs = rangespace([B;Δ])
        f = Fun((x,y)->exp(-x^2-y^2),d)
        @test_throws DimensionMismatch coefficients([0.0;0.0;0.0;0.0;f],rs)
    end

    @testset "off domain evaluate" begin
        g = Fun(1, Segment(Vec(0,-1) , Vec(π,-1)))
        @test g(0.1,-1) ≈ 1
        @test g(0.1,1) ≈ 0

        g = Fun(1, PeriodicSegment(Vec(0,-1) , Vec(π,-1)))
        @test g(0.1,-1) ≈ 1
        @test g(0.1,1) ≈ 0
    end


    @testset "Dirichlet" begin
        testblockbandedoperator(Dirichlet((0..1)^2))
        testblockbandedoperator(Dirichlet((0..1) × (0.0 .. 1)))
        testraggedbelowoperator(Dirichlet(Chebyshev()^2))
        testraggedbelowoperator(Dirichlet(Chebyshev(0..1) * Chebyshev(0.0..1)))
    end
end
