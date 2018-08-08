using ApproxFun, Compat.Test
    import ApproxFun: ChebyshevDirichlet, Ultraspherical, PiecewiseSegment, ContinuousSpace, space,
                        testspace, testbandedoperator, testraggedbelowoperator, testcalculus, testtransforms

@testset "Spaces" begin
    testtransforms(ChebyshevDirichlet{1,1}())

    @test Fun(exp,ChebyshevDirichlet{1,1})(.1) ≈ exp(.1)
    @test Fun(Fun(exp,ChebyshevDirichlet{1,1}),Ultraspherical(1))(.1) ≈ exp(.1)

    d=Interval()
    sp=ChebyshevDirichlet{1,1}(d)
    B=Dirichlet(sp)
    D=Derivative(d)
    L=D^2+I
    @test B[1:2,1:4] ≈ [1 -1 0 0; 1 1 0 0]

    @test csc(2)sin(1 - 0.1)  ≈ ([Dirichlet(d);L]\[[1.,0.],0.])(0.1)
    @test csc(2)sin(1 - 0.1)  ≈ ([B;L]\[[1.,0.],0.])(0.1)

    @test norm(([B;L]\[[1.,0],0])-([Dirichlet(d);L]\[[1.,0],0])) <10eps()




    ## PiecewiseSPace

    x=Fun(identity,Domain(-1..1) \ 0)
    sp=space(x)
    testtransforms(sp;minpoints=2)

    D=Derivative(sp)
    testbandedoperator(D)

    B=[Dirichlet(sp);continuity(sp,0:1)]
    u=[B;
        D^2]\Any[[1,0],zeros(2),0];
    u2=[Dirichlet();Derivative(Chebyshev())^2]\[[1.,0],0]
    @test u(0.) ≈ u2(0.)

    x=Fun(identity,Domain(-10..15) \ [0,1])
    sp=space(x)
    D=Derivative(sp)
    B=Dirichlet(sp)

    u=[B;
        continuity(sp,0:1);
        D^2-x]\[[airyai(-10.),0],zeros(4),0];

    @test u(0.) ≈ airyai(0.)

    s=Fun(sin,-2..2)|>abs
    c=Fun(cos,-2..2)|>abs
    sc=Fun(x->abs(sin(x))+abs(cos(x)),Domain(-2..2) \ [-π/2,0,π/2])
    @test norm(sc-(c+s))<100eps()

    # max/min creates breakpoints

    x=Fun()
    g=4*(x-0.2)
    f=max(-1,g)
    f2=min(f,1)


    @test norm(max(x,x)-x)<100eps()
    @test norm(min(x,x)-x)<100eps()

    f3=Fun(x->x<-0.05?-1.0:(x<0.45?4*(x-.2):1),Domain(-1..1) \ [-0.05,0.45])
    @test norm(f2.(linspace(-1,1,10))-f3.(linspace(-1,1,10))) < 2eps()

    x=Fun(identity, Segment(im,0) ∪ Segment(0,1))
    @test x(0.5) ≈ 0.5
    @test x(0.5im) ≈ 0.5im

    @test Fun(Fun(1.0),space(x))(0.5) == 1.0
    @test Fun(Fun(1.0),space(x))(0.5im) == 1.0

    @test (x+1)(0.5) ≈ 1.5
    @test (x-1)(0.5) ≈ -0.5
    @test (1-x)(0.5) ≈ 0.5



    @test sqrt(1-x)(0.2im) ≈ sqrt(1-0.2im)
    @test sqrt(1-x)(0.2) ≈ sqrt(1-0.2)

    w=2/(sqrt(1-x)*sqrt(1+im*x))
    @test w(0.2im) ≈ 2/(sqrt(1-0.2im)*sqrt(1+im*(0.2im)))
    @test w(0.2) ≈ 2/(sqrt(1-0.2)*sqrt(1+im*(0.2)))

    ## Equivalent spaces

    @test norm(Fun(cos,Chebyshev)-Fun(cos,Jacobi(-0.5,-0.5)))<100eps()
    @test norm(Fun(cos,Chebyshev)-Fun(cos,JacobiWeight(0,0)))<100eps()
    @test norm(Fun(cos,Jacobi(-0.5,-0.5))-Fun(cos,JacobiWeight(0,0))) < 100eps()
    @test norm(Fun(cos,Chebyshev)-Fun(cos,JacobiWeight(0,0,Jacobi(-0.5,-0.5))))<100eps()
    @test norm(Fun(cos,Jacobi(-0.5,-0.5))-Fun(cos,JacobiWeight(0,0,Jacobi(-0.5,-0.5))))<100eps()



    ## ContinuousSpace

    d = PiecewiseSegment(1.,2.,3.,4.)
    S = ContinuousSpace(d)
    testtransforms(S; minpoints=3, invertibletransform=false)

    D = Derivative(S)
    testbandedoperator(D)

    A = [ldirichlet(S); D-I]
    testraggedbelowoperator(A)
    u = [ldirichlet(S); D-I] \ [exp(1.); 0]


    @test u(1.1) ≈ exp(1.1)
    @test u(3.4) ≈ exp(3.4)
    @test last(u) ≈ exp(4)


    d=PiecewiseSegment(0,1.,1.+im,im,0.)
    s=Space(d)

    # narrow down bug
    @test s isa ContinuousSpace
    @test ApproxFun.rangetype(s) == Float64
    cs=ApproxFun.canonicalspace(s)
    @test ApproxFun.rangetype(cs) == Float64

    @test conversion_type(s,cs) == s

    @test Fun(exp,d)(.1) ≈ exp(.1)





    ## Triple SumSpace

    x=Fun()
    w=log(1-x)+sqrt(1-x^2)
    f=w+x
    @test f(0.1) ≈ (w(0.1)+0.1)
    @test (w+1)(0.1) ≈ (w(0.1)+1)
    @test (w+x+1)(0.1) ≈ (w(0.1)+1.1)
    @test ((w+x)+1)(0.1) ≈ (w(0.1)+1.1)


    ## SumSpace bug

    dsp=JacobiWeight(1.,0.,Jacobi(1.,0.,0..1))⊕JacobiWeight(0.5,0.,Jacobi(0.5,-0.5,0..1))
    rsp=Legendre(0..1)⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,0..1))


    C=Conversion(dsp,rsp)

    f=Fun(dsp,[1.,2.,3.,4.,5.])
    @test f(0.1) ≈ (C*f)(0.1)






    ## Piecewise + Cosntant
    Γ=Circle() ∪ Circle(0.0,0.4)
    o=ones(Γ)
    @test o(1.) ≈ 1.0
    @test o(0.4) ≈ 1.0

    G=Fun(z->in(z,component(Γ,2))?[1 0; -1/z 1]:[z 0; 0 1/z],Γ)
    @test (G-I)(exp(0.1im)) ≈ (G(exp(0.1im))-I)


    ## Previoius seffdault

    x=Fun(identity,-1..1)
    f=x+sin(2x)*sqrt(1-x^2)
    @test f(0.1) ≈ 0.1+sin(2*0.1)*sqrt(1-0.1^2)


    ## Check multiple piecewisesapce

    x=Fun(identity,-3 .. -2)+Fun(identity,2..3)
    w=sqrt(9-x^2)
    f=w+Fun()
    @test (f+w)(2.5) ≈ 2w(2.5)
    @test (f+w)(.5) ≈ f(.5)



    ## Check Jacobi recurrence bug

    S=Jacobi(0.5,-0.5)
    f=Fun(exp,S)
    @test f(0.1) ≈ exp(0.1)


    ## Check cancel conversion works
    x=Fun(0..1)
    f=exp(x)-1
    Fun(f,JacobiWeight(1.,0.,0..1))


    ## Arc exp

    z=Fun(identity,Arc(0.,.1,0.,π/2))
    @test exp(z)(0.1exp(0.2im)) ≈ exp(0.1exp(0.2im))



    ## Extending function

    Γ=Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)

    @test isempty(component(Γ,1)\component(Γ,1))
    @test Γ\component(Γ,1) == component(Γ,2) ∪ component(Γ,3)

    @test norm(Fun(ones(component(Γ,1)),Γ) - Fun(x->x ∈ component(Γ,1) ? 1.0 : 0.0,Γ)) == 0


    ## Line

    f=Fun(z->2exp(z^2),PeriodicLine(0.,π/2))
    @test f(1.1im) ≈ 2exp(-1.1^2)


    f=Fun(z->2exp(z^2),Line(0.,π/2))
    @test f(1.1im) ≈ 2exp(-1.1^2)



    ## Exp for Γ

    a=1+10*im;b=2-6*im
    d=Curve(Fun(x->1+a*x+x^2+b*x^3))

    x=Fun(d)

    @test exp(x)(1+a*0.1+0.1^2+b*0.1^3) ≈ exp(1+a*0.1+0.1^2+b*0.1^3)


    ## ChebyshevDirichlet multiplication

    S=ChebyshevDirichlet()
    x=Fun()
    @test norm((ApproxFun.Recurrence(S)*Fun(exp,S)-Fun(x->x*exp(x),S)).coefficients) < 100eps()
    @test norm((x*Fun(exp,S)-Fun(x->x*exp(x),S)).coefficients) < 100eps()


    ## ChebyshevDirichlet Integral
    @test Integral(S,1)*Fun(S,[1.,2.,3.]) ≈ integrate(Fun(Fun(S,[1.,2.,3.]),Chebyshev()))


    ### QuotientSpace test

    import ApproxFun: SpaceOperator

    for (bcs,ret) in ((Dirichlet(Chebyshev()),[1 -1 0 0 0;1 1 0 0 0]),
                      (Neumann(Chebyshev()),[0 1 -4 0 0;0 1 4 0 0]),
                      ([DefiniteIntegral(Chebyshev());SpaceOperator(BasisFunctional(2),Chebyshev(),ConstantSpace())],[2 0 0 0 0;0 1 0 0 0]),
                      (vcat(bvp(Chebyshev(),4)...),[1 -1 1 -1 0;0 1 -4 9 0;1 1 1 1 0;0 1 4 9 0]))
        QS = QuotientSpace(bcs)
        C = Conversion(QS, QS.space)

        norm((bcs*C)[1:size(bcs, 1),1:5] - ret) < 1000eps()
    end



    ## Check union of ChebyshevDirichlet


    dom = Domain(0..1) ∪ Domain(2..3)
    @test components(union(JacobiWeight.(-0.5,-0.5,ChebyshevDirichlet{1,1}.(components(dom)))...)) ==
        (JacobiWeight.(-0.5,-0.5,ChebyshevDirichlet{1,1}.(components(dom)))...)



    ## Ultraspherical special functions

    x = Fun(Ultraspherical(2,0..1))
    sqrt(x)(0.1) ≈ sqrt(0.1)

    f = Fun(x->x*exp(x),Ultraspherical(1,0..1))
    sqrt(f(0.1)) ≈ sqrt(f)(0.1)


    @testset "Fast Ultraspherical*Chebyshev" begin
        f = Fun(exp, Ultraspherical(1))
        g = Fun(exp)


        @test ApproxFun.hasfasttransform(f)
        @test ApproxFun.pointscompatible(f,g)
        @test ApproxFun.hasfasttransformtimes(f,g)
        @test ApproxFun.transformtimes(f,g)(0.1) ≈ exp(0.2)
    end

    @testset "PointSpace" begin
        δ = Fun(ApproxFun.PointSpace([2.0]),[1.0])
        f = Fun(x->cos(50x)) + δ
        g = Fun(x->cos(50x),Ultraspherical(1)) + δ
        @test domain(f) == domain(g)

        @test (f*g)(0.1) ≈ cos(50*0.1)^2
        @test (f*g)(2.0) ≈ 1

        @test exp(f)(0.1) ≈ exp(cos(50*0.1))
        @test exp(f)(2.0) ≈ exp(1)

        @test iszero(sign(-Fun(zero)))
        @test iszero(sign(Fun(zero)))
        @test iszero(abs(-Fun(zero)))
        @test iszero(abs(Fun(zero)))
        @test angle(-Fun(zero)) ≈ Fun(π)
        @test iszero(angle(Fun(zero)))

        @test sign(f)(0.1) ≈ sign(cos(50*0.1))
        @test sign(f)(2.0) ≈ 1
        @test abs(f)(0.1) ≈ abs(cos(50*0.1))
        @test abs(f)(2.0) ≈ 1
        @test angle(f)(0.1) ≈ angle(cos(50*0.1))
        @test angle(f)(2.0) ≈ 0
    end

    x = Fun(Domain(0..1) ∪ Domain(2..3))
    @test length(jumplocations(sign(x))) == 0

    x = Fun(Chebyshev(-1..1))
    @test length(jumplocations(x)) == 0
    @test all(jumplocations(sign(x) + sign(x+0.2)) .≈ [-0.2, 0])

    @testset "ones for SumSpace" begin
        S = Jacobi(0,1) ⊕ JacobiWeight(1/3,0,Jacobi(1/3,2/3)) ⊕ JacobiWeight(2/3,0,Jacobi(2/3,1/3))
        o = ones(S)
        @test o(0.5) == 1
    end


    @testset "blockbandinds for FiniteOperator of pointscompatibleace bug" begin
        S = ApproxFun.PointSpace([1.0,2.0])
        @test ApproxFun.blockbandinds(FiniteOperator([1 2; 3 4],S,S)) == (0,0)
    end
    
    #SumSpace Conversion
    @testset "SumSpace Conversion" begin
        H = HeavisideSpace([-1.0,0.0,1.0])
        C = ContinuousSpace(PiecewiseSegment([-1.0,0,1]))
        S = H + C
        P = Ultraspherical(1,-1.0..0.0) ∪ Ultraspherical(1,0.0..1.0)
        @test Conversion(S,P)[1,1]==1.0
        @test Conversion(S,P)[1,2]==0.5
        @test rangespace(Conversion(S,P))==Ultraspherical(1,-1.0..0.0) ∪ Ultraspherical(1,0.0..1.0)
        @test domainspace(Conversion(S,P))==ApproxFun.SumSpace(HeavisideSpace([-1.0,0.0,1.0]),ContinuousSpace(PiecewiseSegment([-1.0,0,1])))
    end
end
