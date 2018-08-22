using ApproxFun, SpecialFunctions, Random, Test
    import ApproxFun: HeavisideSpace, PointSpace, DiracSpace, PiecewiseSegment

@testset "Singularities" begin
    @testset "sqrt" begin
        x=Fun(identity);
        @test sqrt(cos(π/2*x))(.1) ≈ sqrt(cos(.1π/2))

        x=Fun(identity,-2..2)
        u=sqrt(4-x^2)/(2π)

        @test u(.1) ≈ sqrt(4-0.1^2)/(2π)
        @test sum(u) ≈ 1

        #this call threw an error, which we check
        @test length(values(u)) == 1


        f = Fun(x->x*cot(π*x/2))
        x = Fun(identity)
        u = Fun(JacobiWeight(1.,1.,Interval()), (f/(1-x^2)).coefficients)
        @test 1/(0.1*cot(π*.1/2)) ≈ (1/u)(.1)

        @test (x/u)(.1) ≈ tan(π*.1/2)

        f = Fun(x->exp(-x^2),Line(0.,0.,-.5,-.5),400)
        @test sum(f) ≈ sqrt(π)

        f=Fun(x->exp(x)/sqrt(1-x.^2),JacobiWeight(-.5,-.5))
        @test f(.1) ≈ (x->exp(x)/sqrt(1-x.^2))(.1)

        @test norm(Fun(exp,Legendre(0..1))+sqrt(Fun(0..1))) ≈ 2.491141949903508

        # sampling Chebyshev
        x=Fun(identity)
        f = exp(x)/sqrt(1-x^2)
        g = cumsum(f)
        @test abs(g(-1)) ≤ 1E-15
        @test g'(0.1) ≈ f(0.1)
    end

    @testset "JacobiWeight Derivative" begin
        S=JacobiWeight(-1.,-1.,Chebyshev(0..1))

        # Checks bug in Derivative(S)
        @test typeof(ConstantSpace(Domain(0..1))) <: Space{Segment{Float64},Float64}

        D=Derivative(S)
        f=Fun(S,Fun(exp,0..1).coefficients)
        x=0.1
        @test f(x) ≈ exp(x)*x^(-1)*(1-x)^(-1)/4
        @test (D*f)(x) ≈ -exp(x)*(1+(x-3)*x)/(4*(x-1)^2*x^2)


        S=JacobiWeight(-1.,0.,Chebyshev(0..1))
        D=Derivative(S)

        f=Fun(S,Fun(exp,0..1).coefficients)
        x=.1
        @test f(x) ≈ exp(x)*x^(-1)/2
        @test (D*f)(x) ≈ exp(x)*(x-1)/(2x^2)
    end


    ## ODEs

    ## f/g bugs

    @testset "Jacobi singularity" begin
        x = Fun(identity)
        f = exp(x)/(1-x.^2)

        @test f(.1) ≈ exp(.1)/(1-.1^2)
        f = exp(x)/(1-x.^2).^1
        @test f(.1) ≈ exp(.1)/(1-.1^2)
        f = exp(x)/(1-x.^2).^1.0
        @test f(.1) ≈ exp(.1)/(1-.1^2)



        ## 1/f with poles

        x=Fun(identity)
        f=sin(10x)
        g=1/f

        @test g(.123) ≈ csc(10*.123)
    end

    @testset "Ray and Line" begin
        @test Inf in Ray()   # this was a bug

        f=Fun(x->exp(-x),0..Inf)
        @test f'(.1) ≈ -f(.1)

        x=Fun(identity,Ray())
        f=exp(-x)
        u=integrate(f)
        @test (u(1.)-u(0)-1) ≈ -f(1)

        x=Fun(identity,Ray())
        f=x^(-0.123)*exp(-x)
        @test integrate(f)'(1.) ≈ f(1.)

        @test ≈(sum(Fun(sech,0..Inf)),sum(Fun(sech,0..40));atol=1000000eps())

        f=Fun(sech,Line())
        Fun(f,Ray())(2.0) ≈ sech(2.0)
        Fun(f,Ray(0.,π))(-2.0) ≈ sech(-2.0)
        Fun(sech,Ray(0.,π))(-2.0) ≈ sech(-2.0)


        #Ei (Exp Integral)

        y=Fun(Ray())
        q=integrate(exp(-y)/y)
        @test (q-last(q))(2.) ≈ (-0.04890051070806113)

        ## Line
        f=Fun(x->exp(-x^2),Line())

        @test f'(0.1) ≈ -2*0.1exp(-0.1^2)
        @test (Derivative()*f)(0.1) ≈ -2*0.1exp(-0.1^2)

        ## PeriodicLine

        d=PeriodicLine()
        D=Derivative(d)

        f = Fun(x->sech(x-0.1),d,200)
        @test f(1.) ≈ sech(1-0.1)


        f=Fun(x->sech(x-0.1),d)
        @test f(1.) ≈ sech(1-0.1)

        @test ≈((D*f)(.2),-0.0991717226583897;atol=100000eps())
        @test ≈((D^2*f)(.2),-0.9752522555114987;atol=1000000eps())
    end

    @testset "LogWeight" begin
        x=Fun(identity,-1..1)
        f=exp(x+1)-1
        @test log(f)(0.1) ≈ log(f(0.1))

        x=Fun(identity,0..1)
        f=exp(x)-1
        @test log(f)(0.1) ≈ log(f(0.1))

        x=Fun(identity,0..1)
        @test Fun(exp(x)/x-1/x,Chebyshev)(0.1) ≈ (exp(0.1)-1)/0.1

        x=Fun(identity,0..1)
        f=1/x
        p=integrate(f)
        @test (p-p(1.))(0.5) ≈ log(0.5)

        f=1/(1-x)
        p=integrate(f)
        @test (p-p(0.))(0.5) ≈ -log(1-0.5)
    end

    @testset "Complex domains sqrt" begin
        a=1+10*im;b=2-6*im
        d=Curve(Fun(x->1+a*x+b*x^2))

        x=Fun(d)
        w=sqrt(abs(first(d)-x))*sqrt(abs(last(d)-x))

        @test sum(w/(x-2.))/(2π*im) ≈ (-4.722196879007759+2.347910413861846im)
        @test linesum(w*log(abs(x-2.)))/π ≈ (88.5579588360686)

        a=Arc(0.,1.,0.,π/2)
        ζ=Fun(identity,a)
        f=Fun(exp,a)*sqrt(abs((ζ-1)*(ζ-im)))
    end

    @testset "DiracDelta and PointSpace" begin
        a,b=DiracDelta(0.),DiracDelta(1.)
        f=Fun(exp)
        g=a+0.2b+f
        @test components(g)[2](0.) ≈ 1.
        @test g(.1) ≈ exp(.1)
        @test sum(g) ≈ (sum(f)+1.2)

        #Checks prevoius bug
        δ=DiracDelta()
        x=Fun()
        w=sqrt(1-x^2)
        w+δ


        ## PointSpace

        @test eltype(domain(PointSpace([0,0.1,1])) ) == Float64

        f=Fun(x->(x-0.1),PointSpace([0,0.1,1]))
        @test roots(f) == [0.1]

        a=Fun(exp,space(f))
        @test f/a == Fun(x->(x-0.1)*exp(-x),space(f))

        f = Fun(space(f),[1.,2.,3.])

        g = f + Fun(2..3)
        @test f(0.0) ≈ g(0.0) ≈ 1.0
        @test f(0.1) ≈ g(0.1) ≈ 2.0
        @test f(1.0) ≈ g(1.0) ≈ 3.0

        @test g(2.3) ≈ 2.3


        h = a + Fun(2..3)

        # for some reason this test is broken only on Travis
        @test_skip g/h ≈ f/a + Fun(1,2..3)
    end

    @testset "DiracDelta integration and differentiation" begin
        δ = DiracDelta()
        h = integrate(δ)
        @test domain(h) == PiecewiseSegment([0,Inf])
        @test h(-2) == 0
        @test h(2) == 1

        δ = 0.3DiracDelta(0.1) + 3DiracDelta(2.3)
        h = integrate(δ)
        @test domain(h) == PiecewiseSegment([0.1,2.3,Inf])
        @test h(-2) == 0
        @test h(2) == 0.3
        @test h(3) == 3.3

        δ = (0.3+1im)DiracDelta(0.1) + 3DiracDelta(2.3)
        h = integrate(δ)
        @test domain(h) == PiecewiseSegment([0.1,2.3,Inf])
        @test h(-2) == 0
        @test h(2) == 0.3+1im
        @test h(3) == 3.3+1im
    end

    @testset "DiracDelta sampling" begin
        δ = 0.3DiracDelta(0.1) + 3DiracDelta(2.3)
        Random.seed!(0)
        for _=1:10
            @test sample(δ) ∈ [0.1, 2.3]
        end
        Random.seed!(0)
        r = sample(δ, 10_000)
        @test count(i -> i == 0.1, r)/length(r) ≈ 0.3/(3.3) atol=0.01
    end

    @testset "Multiple roots" begin
        x=Fun(identity,-1..1)
        @test (1/x^2)(0.1) ≈ 100.
        @test (1/x^2)(-0.1) ≈ 100.

        fc=x*(1+x)^2
        @test (1/fc)(0.1) ≈ 1/fc(0.1)

        fc=x*(1-x)^2
        @test (1/fc)(0.1) ≈ 1/fc(0.1)
    end

    @testset "special function singularities" begin
        x=Fun(0..1)
        @test erf(sqrt(x))(0.1) ≈ erf(sqrt(0.1))
        @test erfc(sqrt(x))(0.1) ≈ erfc(sqrt(0.1))

        ## roots of log(abs(x-y))
        x=Fun(-2..(-1))
        @test roots(abs(x+1.2)) ≈ [-1.2]

        f=abs(x+1.2)

        @test norm(abs(f)-f)<10eps()
        @test norm(sign(f)-Fun(1,space(f)))<10eps()


        @test log(f)(-1.3) ≈ log(abs(-1.3+1.2))
        @test log(f)(-1.1) ≈ log(abs(-1.1+1.2))

        #393
        x=Fun(0..1)
        f = exp(x)*sqrt(x)*log(1-x)
        @test f(0.1) ≈ exp(0.1)*sqrt(0.1)*log(1-0.1)
    end


    @testset "Jacobi conversions" begin
        S1,S2=JacobiWeight(3.,1.,Jacobi(1.,1.)),JacobiWeight(1.,1.,Jacobi(0.,1.))
        f=Fun(S1,[1,2,3.])
        C=Conversion(S1,S2)
        Cf=C*f
        @test Cf(0.1) ≈ f(0.1)

        S1,S2=JacobiWeight(3.,2.,Jacobi(1.,1.)),JacobiWeight(1.,1.,Jacobi(0.,0.))
        f=Fun(S1,[1,2,3.])
        C=Conversion(S1,S2)
        Cf=C*f
        @test Cf(0.1) ≈ f(0.1)
    end


    @testset "Derivative operator for HeavisideSpace" begin
        H = HeavisideSpace([-1.0,0.0,1.0])
        @test Fun(H, [1.0])(1.0) == 0.0
        @test Fun(H, [0.0,1.0])(1.0) == 1.0

        H=HeavisideSpace([1,2,3])
        D=Derivative(H)
        @test domain(D)==PiecewiseSegment([1,2,3])
        @test D[1,1]==-1
        @test D[1,2]==1

        H=HeavisideSpace([1,2,3,Inf])
        D=Derivative(H)
        @test domain(D)==PiecewiseSegment([1,2,3,Inf])
        @test D[1,1]==-1
        @test D[2,2]==-1
        @test D[1,2]==1

        S = HeavisideSpace([-1.0,0.0,1.0])
        @test Derivative(S) === Derivative(S,1)
    end
end
