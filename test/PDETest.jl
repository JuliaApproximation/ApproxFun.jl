using ApproxFun, LinearAlgebra, Test
    import ApproxFun: testbandedblockbandedoperator, testblockbandedoperator, testraggedbelowoperator, Block

@testset "PDE" begin
    @testset "Zero Dirichlet" begin
        S = JacobiWeight(1.,1.,Jacobi(1.,1.))^2
        Δ = Laplacian(S)

        testbandedblockbandedoperator(Δ)

        u = Fun((x,y)->sin(π*x)*sin(π*y),S)
        f = -2π^2*u

        F = qr(Δ)
        ApproxFun.resizedata!(F,:,1000)
        @time v=F\f
        @test norm((u-v).coefficients)<100eps()


        F=qr(Δ)
        ApproxFun.resizedata!(F.R_cache,:,100)
        ApproxFun.resizedata!(F.R_cache,:,1000)
        @time v=F \ f
        @test norm((u-v).coefficients)<100eps()

        F=qr(Δ)
        @time v=F\f
        @test norm((u-v).coefficients)<100eps()
    end

    @testset "Rectangle Laplace/Poisson" begin
        dx = dy = ChebyshevInterval()
        d = dx × dy
        g = Fun((x,y)->exp(x)*cos(y),∂(d))

        B = Dirichlet(d)

        testblockbandedoperator(B)
        testbandedblockbandedoperator(Laplacian(d))
        testbandedblockbandedoperator(Laplacian(d)+0.0I)

        A=[Dirichlet(d);Laplacian(d)]

        @time u=A\[g,0.]
        @test u(.1,.2) ≈ real(exp(0.1+0.2im))

        A=[Dirichlet(d);Laplacian(d)+0.0I]
        @time u=A\[g,0.]

        @test u(.1,.2) ≈ real(exp(0.1+0.2im))

        ## Poisson
        f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),ChebyshevInterval()^2,500)  #default is [-1,1]^2
        d=domain(f)
        A=[Dirichlet(d);Laplacian(d)]
        @time  u=\(A,[zeros(∂(d));f];tolerance=1E-7)
        @test ≈(u(.1,.2),-0.04251891975068446;atol=1E-5)
    end

    @testset "Bilaplacian" begin
        dx = dy = ChebyshevInterval()
        d = dx × dy
        Dx = Derivative(dx); Dy = Derivative(dy)
        L = Dx^4⊗I + 2*Dx^2⊗Dy^2 + I⊗Dy^4

        testbandedblockbandedoperator(L)

        B = Dirichlet(dx) ⊗ Operator(I,dy)
        testraggedbelowoperator(B)

        A=[Dirichlet(dx) ⊗ Operator(I,dy);
                Operator(I,dx)  ⊗ Dirichlet(dy);
                Neumann(dx) ⊗ Operator(I,dy);
                Operator(I,dx) ⊗ Neumann(dy);
                 L]

        testraggedbelowoperator(A)

        @time u=\(A,[[1,1],[1,1],[0,0],[0,0],0];tolerance=1E-5)
        @test u(0.1,0.2) ≈ 1.0
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

    @testset "Schrodinger" begin
        dx=0..1; dt=0.0..0.001
        C=Conversion(Chebyshev(dx)*Ultraspherical(1,dt),Ultraspherical(2,dx)*Ultraspherical(1,dt))
        testbandedblockbandedoperator(C)
        testbandedblockbandedoperator(Operator{ComplexF64}(C))

        d = dx × dt

        x,y = Fun(d)
        @test x(0.1,0.0001) ≈ 0.1
        @test y(0.1,0.0001) ≈ 0.0001

        V = x^2

        Dt=Derivative(d,[0,1]);Dx=Derivative(d,[1,0])

        ϵ = 1.
        u0 = Fun(x->exp(-100*(x-.5)^2)*exp(-1/(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
        L = ϵ*Dt+(.5im*ϵ^2*Dx^2)
        testbandedblockbandedoperator(L)

        @time u = \([timedirichlet(d);L],[u0,[0.,0.],0.];tolerance=1E-5)
        @test u(0.5,0.001) ≈ 0.857215539785593+0.08694948835021317im  # empircal from ≈ schurfact
    end

    @testset "check we dispatch correctly to get fast build" begin
        S = JacobiWeight(1.,1.,Jacobi(1.,1.))^2
        Δ = Laplacian(S)
        @time S = view(Δ.op.ops[1].ops[1].op,Block.(1:40), Block.(1:40))
        @test typeof(S.parent.domaintensorizer) == ApproxFun.Trivial2DTensorizer
    end
end
