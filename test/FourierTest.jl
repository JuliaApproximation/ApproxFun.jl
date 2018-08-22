using ApproxFun, Test, BlockArrays, BlockBandedMatrices, SpecialFunctions
    import ApproxFun: testspace, testtransforms, testmultiplication,
                      testbandedoperator, testcalculus, Block
    import SpecialFunctions: factorial



@testset "Fourier" begin
    @testset "Fourier transforms" begin
        for d in (PeriodicInterval(0.1,0.5),Circle(1.0+im,2.0))
            testtransforms(CosSpace(d);minpoints=2)
            testtransforms(SinSpace(d))

            testtransforms(Taylor(d))
            testtransforms(Hardy{false}(d))

            testspace(Laurent(d);hasintegral=false)
            testspace(Fourier(d);hasintegral=false)
        end

        P = ApproxFun.plan_transform(Fourier(),4)
        v = randn(4) .+ im.*randn(4)
        @test P*v == P*real(v) + im*(P*imag(v))
    end


    @test sum(Fun(1,CosSpace())) ≈ 2π
    @test sum(Fun(SinSpace(),[1])) == 0


    f=Fun(t->cos(t)+cos(3t),CosSpace)
    @test f(0.1) ≈ cos(0.1)+cos(3*0.1)
    @test (f*f-Fun(t->(cos(t)+cos(3t))^2,CosSpace)).coefficients|>norm <100eps()


    f=Fun(exp,Taylor(Circle()))
    @test f(exp(0.1im)) ≈ exp(exp(0.1im))
    @test f(1.0) ≈ exp(1.0)
    g=Fun(z->1/(z-0.1),Hardy{false}(Circle()))
    @test (f(1.)+g(1.)) ≈ (exp(1.) + 1/(1-.1))


    ## Periodic
    f=Fun(x->exp(-10sin((x-.1)/2)^2),Laurent)
    @test f(.5) ≈ (Conversion(space(f),Fourier(domain(f)))*f)(.5)
    @test f(.5) ≈ Fun(f,Fourier)(.5)


    Γ=Circle(1.1,2.2)
    z=Fun(Fourier(Γ))
    @test space(z)==Fourier(Γ)
    @test z(1.1+2.2exp(0.1im)) ≈ 1.1+2.2exp(0.1im)




    @test Fun(Laurent(0..2π),[1,1.,1.])(0.1) ≈ 1+2cos(0.1)
    @test Fun(Laurent(-1..1),[1,1.,1.])(0.1) ≈ 1+2cos(π*(0.1+1))
    @test Fun(Laurent(0..1),[1,1.,1.])(0.1) ≈ 1+2cos(2π*0.1)


    @test abs(Fun(cos,Circle())(exp(0.1im))-cos(exp(0.1im)))<100eps()
    @test abs(Fun(cos,Circle())'(exp(0.1im))+sin(exp(0.1im)))<100eps()
    @test abs(Fun(cos,Circle())'(exp(0.1im))+Fun(sin,Circle())(exp(0.1im)))<100eps()

    @test norm(Fun(x->Fun(cos,Fourier(-π .. π),20)(x),20)-Fun(cos,20)) <100eps()
    @test norm(Fun(x->Fun(cos,Fourier(-π .. π))(x))-Fun(cos)) <100eps()
    @test norm(Fun(cos,Fourier)'+Fun(sin,Fourier)) < 100eps()
    @test norm(Fun(x->Fun(cos,Laurent)(x))-Fun(cos)) <100eps()
    @test norm(Fun(cos,Laurent)'+Fun(sin,Laurent)) < 100eps()
    @test norm(Fun(cos,Circle())'+Fun(sin,Circle()))<100eps()



    for f in (Fun(θ->sin(sin(θ)),SinSpace()),Fun(θ->cos(θ)+cos(3θ),CosSpace()),
                Fun(θ->sin(sin(θ)),Fourier()),Fun(θ->cos(θ)+cos(3θ),CosSpace()))
        @test norm(integrate(f)'-f)<10eps()
    end

    let f=Fun(exp,Circle())
        @test component(f,1)(exp(0.1im)) ≈ exp(exp(0.1im))
        @test f(exp(0.1im)) ≈ exp(exp(0.1im))
        @test norm(f'-f)<100eps()
        @test norm(integrate(f)+1-f)<100eps()
    end

    let f=Fun(x->exp(-10sin((x-.1)/2)^2),Fourier)
        @test real(f)(.1) ≈ f(.1)
    end



    @test (Fun(z->sin(z)*cos(1/z),Circle())*Fun(z->exp(z)*airyai(1/z),Circle()))(exp(.1im)) ≈
                    (z->sin(z)*cos(1/z)*exp(z)*airyai(1/z))(exp(.1im))

    ## Calculus

    let f=Fun(t->cos(t),CosSpace)
        D=Derivative(space(f))
        @test (D*f)(.1) ≈ -sin(.1)
        @test f'(.1) ≈ -sin(.1)
    end

    let f=Fun(t->sin(t),SinSpace)
        D=Derivative(space(f))
        @test (D*f)(.1) ≈ cos(.1)
        @test f'(.1) ≈ cos(.1)
    end

    let f=Fun(cos,Fourier)
        @test norm((Derivative(space(f))^2)*f+f)<100eps()
    end



    ## Taylor


    @test Fun(Taylor())  == Fun(Taylor(),[0.,1.])

    @test Fun(Taylor())(1.0) ≈ 1.0
    @test Fun(Taylor(Circle(0.1,2.2)))(1.0) ≈ 1.0
    @test Fun(Taylor(Circle(0.1+0.1im,2.2)))(1.0) ≈ 1.0


    @test Multiplication(Fun(Taylor()),Taylor())[1:3,1:3] == [0. 0. 0.; 1. 0. 0.; 0. 1. 0.]

    for d in (Circle(),Circle(0.5),Circle(-0.1,2.))
        S=Taylor(d)
        D=Derivative(S)
        ef=Fun(exp,S)
        @test norm((D*ef-ef).coefficients)<4000eps()
        @test norm((D^2*ef-ef).coefficients)<200000eps()
        u=[Evaluation(S,0.),D-I]\[1.;0.]
        @test norm((u-ef).coefficients)<200eps()
        @test norm((Integral(S)*Fun(exp,S)+ef.coefficients[1]-ef).coefficients)<100eps()


        f=Fun(z->exp(1/z)-1,Hardy{false}(d))
        df=Fun(z->-1/z^2*exp(1/z),Hardy{false}(d))
        @test norm((Derivative()*f-df).coefficients)<1000eps()
        @test norm((Derivative()^2*f-df').coefficients)<100000eps()
        @test norm((f'-df).coefficients)<1000eps()
    end

    d=Circle()
    S=Taylor(d)
    D=Derivative(S)
    D-I
    ef=Fun(exp,S)
    @test norm((D*ef-ef).coefficients)<1000eps()
    @test norm((D^2*ef-ef).coefficients)<100000eps()
    u=[Evaluation(S,0.),D-I]\[1.;0.]

    # check's Derivative constructor works
    D=Derivative(Taylor(PeriodicInterval()))





    ## Multiplication

    s=Fun(t->(sin(t)+sin(2t))*cos(sin(t)),SinSpace)
    b=Fun(t->(sin(t)+sin(3t)),SinSpace)

    @test (s*s)(.1) ≈ s(.1)^2
    @test (s*b)(.1) ≈ s(.1)*b(.1)

    s=Fun(t->(cos(t)+cos(2t))*cos(cos(t)),CosSpace)
    b=Fun(t->(1+cos(t)+cos(3t)),CosSpace)

    @test (s*s)(.1) ≈ s(.1)^2
    @test (s*b)(.1) ≈ s(.1)*b(.1)

    s=Fun(t->(cos(t)+cos(2t))*cos(cos(t)),CosSpace)
    b=Fun(t->(sin(t)+sin(3t)),SinSpace)

    @test (s*b)(.1) ≈ s(.1)*b(.1)


    s=Fun(t->(sin(t)+sin(2t))*cos(sin(t)),SinSpace)
    b=Fun(t->(1+cos(t)+cos(3t)),CosSpace)

    @test (s*b)(.1) ≈ s(.1)*b(.1)



    a=Fun(t->exp(cos(t)+sin(t)),Fourier)
    b=Fun(t->sin(t)+cos(3t)+1,Fourier)

    @test (a*b)(.1) ≈ a(.1)*b(.1)

    a=Fun(t->exp(cos(t)),CosSpace)
    b=Fun(t->sin(t)+cos(3t)+1,Fourier)

    @test (a*b)(.1) ≈ a(.1)*b(.1)

    a=Fun(t->sin(sin(t)),SinSpace)
    b=Fun(t->sin(t)+cos(3t)+1,Fourier)

    @test (a*b)(.1) ≈ a(.1)*b(.1)




    # Check bug in off centre Circle
    c2=-0.1+.2im;r2=0.3;
    d2=Circle(c2,r2)
    z=Fun(identity,d2)

    @test z(-0.1+.2im+0.3*exp(0.1im)) ≈ (-0.1+.2im+0.3*exp(0.1im))



    # false Circle
    @test Fun(exp,Fourier(Circle(0.,1.,false)))(exp(0.1im)) ≈ exp(exp(.1im))
    @test Fun(exp,Laurent(Circle(0.,1.,false)))(exp(0.1im)) ≈ exp(exp(.1im))



    ## Reverse orientation

    f=Fun(z->1/z,Taylor(1/Circle()))
    @test f(exp(0.1im)) ≈ exp(-0.1im)



    ## exp(z)

    z=Fun(identity,Circle())
    cfs=exp(z).coefficients[1:2:end]
    for k=1:length(cfs)
        @test abs(cfs[k]-1/factorial(1.0(k-1))) ≤ 1E-10
    end



    ##  Norms


    @test sum(Fun(CosSpace(),[1.]))/(2π) ≈ 1.
    @test sum(Fun(CosSpace(),[0.,1.])^2)/(2π) ≈ 0.5
    @test sum(Fun(CosSpace(),[0.,0.,1.])^2)/(2π) ≈ 0.5
    @test sum(Fun(CosSpace(),[0.,0.,0.,1.])^2)/(2π) ≈ 0.5


    @test sum(Fun(SinSpace(),[0.,1.])^2)/(2π) ≈ 0.5
    @test sum(Fun(SinSpace(),[0.,0.,1.])^2)/(2π) ≈ 0.5
    @test sum(Fun(SinSpace(),[0.,0.,0.,1.])^2)/(2π) ≈ 0.5


    ## Bug in multiplicaiton

    @test Fun(SinSpace(),Float64[])^2 == Fun(SinSpace(),Float64[])
    @test Fun(Fourier(),[1.])^2 ≈ Fun(Fourier(),[1.])

    B=Evaluation(Laurent(0..2π),0,1)
    @test B*Fun(sin,domainspace(B)) ≈ 1.0



    ## Conversion between reverse

    C = Conversion(SinSpace()⊕CosSpace(),Fourier())
    @test C[Block(1), Block(1)] ≈ [0 1; 1 0]
    @test ApproxFun.defaultgetindex(C, Block.(1:2), Block.(1:2)) isa AbstractMatrix

    testbandedoperator(C)

    ## Diagonal Derivative
    D = Derivative(Laurent())
    @test isdiag(D)


    ## Test bug in multiplication

    y = Fun(Circle())
    @test (y^2) ≈ Fun(z->z^2,domain(y))


    ## Test Multiplication
    mySin = Fun(Fourier(),[0,1.0])
    A = Multiplication(mySin,Fourier())
    @test A.op[1,1] == 0


    mySin = Fun(Fourier(),[0,1])
    A = Multiplication(mySin,Fourier())
    @test A.op[1,1] == 0



    @testset "Negatively oriented circles" begin
        f1 = Fun(z -> exp(1/z), Circle(0.0,0.2))
        f̃1 = Fun(z -> exp(1/z), Circle(0.0,0.2,false))
        f̃2 = Fun(z -> exp(1/z), Circle(0.0,0.3,false))

        @test f1(0.2exp(0.1im)) ≈ exp(1/(0.2exp(0.1im)))
        @test f̃1(0.2exp(0.1im)) ≈ exp(1/(0.2exp(0.1im)))
        @test f̃2(0.3exp(0.1im)) ≈ exp(1/(0.3exp(0.1im)))

        @test sum(f1) ≈ -sum(f̃1)
        @test sum(f̃1) ≈ sum(f̃2)
    end


    @testset "Fourier inplace" begin
        S = Fourier()

        x = [1.,2,3,4,5]
        y = similar(x)
        z = similar(x)
        P = ApproxFun.plan_transform(S, x)
        P! = ApproxFun.plan_transform!(S, x)
        mul!(y, P, x)
        @test x ≈ [1.,2,3,4,5]
        mul!(z, P!, x)
        @test x ≈ [1.,2,3,4,5]
        @test y ≈ z ≈ P*x ≈ P!*copy(x)

        P = ApproxFun.plan_itransform(S, x)
        P! = ApproxFun.plan_itransform!(S, x)
        mul!(y, P, x)
        @test x ≈ [1.,2,3,4,5]
        mul!(z, P!, x)
        @test x ≈ [1.,2,3,4,5]
        @test y ≈ z ≈ P*x ≈ P!*copy(x)
    end
end
