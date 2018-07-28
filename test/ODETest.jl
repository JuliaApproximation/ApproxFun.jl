using ApproxFun, Domains, SpecialFunctions, Test
    import ApproxFun: Multiplication, testraggedbelowoperator, testbandedoperator, interlace, ∞

@testset "ODE" begin
    @testset "Airy" begin
        d=Interval(-10.,5.);
        S=Chebyshev(d)


        Bm=Evaluation(d,d.a);
        Bp=Evaluation(d,d.b);
        B=[Bm;Bp];
        D2=Derivative(d,2);
        X=Multiplication(Fun(x->x,d));

        testbandedoperator(D2-X)
        testraggedbelowoperator([B;D2-X])

        @time u = [B;D2-X] \ [airyai(d.a),airyai(d.b),0.];
        @test Number.(Array(B*u)) ≈ [airyai(d.a),airyai(d.b)]

        @test ≈(u(0.),airyai(0.);atol=10ncoefficients(u)*eps())

        @time u=[Bm;D2-X;Bp]\[airyai(d.a),0.,airyai(d.b)];
        @test ≈(u(0.),airyai(0.);atol=10ncoefficients(u)*eps())

        @time u=[D2-X;Bm;Bp]\[0.,airyai(d.a),airyai(d.b)];
        @test ≈(u(0.),airyai(0.);atol=10ncoefficients(u)*eps())

        d=Interval(-1000.,5.);
        Bm=Evaluation(d,d.a);
        Bp=Evaluation(d,d.b);
        B=[Bm;Bp];
        D2=Derivative(d,2);
        X=Multiplication(Fun(x->x,d));

        u=[B;D2-X]\[airyai(d.a),airyai(d.b),0.];
        @test ≈(u(0.),airyai(0.);atol=10ncoefficients(u)*eps())

        B=Neumann(d);
        A=[B;D2-X];
        b=[[airyaiprime(d.a),airyaiprime(d.b)],0.];

        @time u=A\b;

        @test ≈(u(0.),airyai(0.);atol=10ncoefficients(u)*eps())

        ##) Neumann condition
    end

    @testset "Exp" begin
        f=Fun(x->-x^2)
        g=Fun(t->exp(-t^2))

        fp=f';
        Bm=Evaluation(domain(f),leftendpoint(domain(f)));
        u=[Bm,Derivative(domain(f)) - fp]\[exp(f(leftendpoint(domain(f)))),0.];
        @test norm(u-g)<100eps()
    end

    @testset "Oscillatory integral" begin
        f=Fun(exp);
        D=Derivative(domain(f));
        w=10.;
        B=ApproxFun.SpaceOperator(BasisFunctional(floor(w)),Chebyshev(),ApproxFun.ConstantSpace(Float64));
        A=[B; D+1im*w*I];

        @time u = A\[0.,f];
        @test (u(1.)exp(1im*w)-u(-1.)exp(-1im*w)) ≈ (-0.18575766879136255 + 0.17863980562549928im)
    end

    @testset "Bessel" begin
        d=ChebyshevInterval()
        D=Derivative(d)
        x=Fun(identity,d)
        A=x^2*D^2+x*D+x^2
        testbandedoperator(x^2*D^2)
        testbandedoperator(ToeplitzOperator([0.5],[0.0,0.5]))
        testbandedoperator(HankelOperator(Float64[]))
        testbandedoperator(A)
        u=[ldirichlet(d);A]\[besselj(0,leftendpoint(d)),0.];

        @test u(0.1) ≈ besselj(0.,0.1)
        @test norm(A*u)<10eps()
        @test norm(Fun(A.ops[1]*u,d)-x.^2 .* differentiate(u,2))<eps()
        @test norm(Fun(A.ops[2]*u,d)-x.*u') < eps()
        @test norm(Fun(A.ops[end]*u,d)-x.^2 .* u) < eps()
        @test norm(x.^2 .*u'' + x.*u' + x.^2 .* u)<10eps()
    end

    @testset "QR" begin
        S=Chebyshev()
        B=Dirichlet(S)
        D=Derivative(S)

        Q,R=qr([B;D^2+I])
        @test Q[1,1] == -0.5773502691896257
        @test size(Q') == (∞,∞)
        u=R\(Q'*[[cos(-1.0),cos(1.0)],0.0])

        @test u(0.) ≈ cos(0.0)

        S=Chebyshev()
        A=[Dirichlet(S);Derivative(S)^2 - I]
        QRf = qr(A)
        @test (QRf\[[1.,0],0])(0.0) ≈ 0.3240271368319427
        Q,R = qr(A)
        u=(R\(Q'*[[1.,0.0],0.0]))
        @test u(0.0)  ≈ 0.3240271368319427
        # matrix RHS
        U = QRf \ [[[1.,0.],0] [[0.,1.0],0.]]
        @test U[1,1] ≈ u

        x=Fun(S)
        A=[Dirichlet(S); Derivative(S)^2 - exp(im*x)]
        QRf = qr(A)

        u=(QRf\[[1.,0.0],0.0])
        @test u(0.0) ≈ (0.3329522068795961 + 0.024616008954634165im)
    end

    @testset "Union of intervals" begin
        x=Fun(identity,Domain(-2..15) \ [-1,0])
        sp=space(x)

        B = [Dirichlet(sp);continuity(sp,0:1)]

        # We don't want to concat piecewise space
        @test !(continuity(sp,0) isa ApproxFun.VectorInterlaceOperator)
        @test B isa ApproxFun.VectorInterlaceOperator

        D=Derivative(sp)
        A=[B;D^2-x]

        ApproxFun.testraggedbelowoperator(A)
        QR=qr(A)

        @time u=QR\[[airyai(-2.),0.0],zeros(4),0.0]

        @test u(0.0) ≈ airyai(0.)
    end

    @testset "Vector" begin
        d=ChebyshevInterval()
        D=Derivative(d);
        B=ldirichlet();
        Bn=lneumann();

        f=Fun(x->[exp(x),cos(x)],d)

        A=[B 0;
           Bn 0;
           0 B;
           D^2-I 2.0I;
           0 D+I]

        # makes sure ops are in right order
        @test A.ops[4,1] isa ApproxFun.PlusOperator
        QR=qr(A)
        v=Any[0.,0.,0.,f...]
        @test (QR\v)(0.0) ≈ [0.0826967758420519,0.5553968826533497]


        Q,R=qr(A)
        v=Any[0.,0.,0.,f...]
        @test (QR\v)(0.0) ≈ [0.0826967758420519,0.5553968826533497]
    end

    @testset "Auto-space" begin
        t=Fun(identity,0..1000)
        L=𝒟^2+2I  # our differential operator, 𝒟 is equivalent to Derivative()

        u=[ivp();L]\[0.;0.;cos(100t)]
        @test ≈(u(1000.0),0.00018788162639452911;atol=1000eps())

        x=Fun(identity,1..2000)
        d=domain(x)
        B=Dirichlet()
        ν=100.
        L=(x^2*𝒟^2) + x*𝒟 + (x^2 - ν^2)   # our differential operator

        @time u=[B;L]\[[besselj(ν,first(d)),besselj(ν,last(d))],0.]
        @test ≈(u(1900.),besselj(ν,1900.);atol=1000eps())
    end

    @testset "Complex RHS for real operator" begin
        D=Derivative(Chebyshev())
        B=ldirichlet()
        u1=[B;D]\[0.;Fun(exp)+0im]
        u2=[B;D]\[0.;Fun(exp)]
        @test u1(0.1) ≈ u2(0.1)
    end

    @testset "Delay" begin
        S = Chebyshev(0..1)
        D = Derivative(S)
        L = [D            0I      0I;
             I             D      0I;
             0I           I        D]


        B = [ldirichlet(S)   0                 0;
             rdirichlet(S)  -ldirichlet(S)     0;
             0              rdirichlet(S)      -ldirichlet(S)]

        v = [B; L] \ [2.3; zeros(2); zeros(3)]

        u = v[1] + setdomain(v[2], Domain(1..2)) + setdomain(v[3], Domain(2..3))

        @test abs(u'(1.3) + u(1.3-1) ) ≤ 10eps()
    end
end
