using ApproxFun, Base.Test

@test norm(Fun(x->Fun(cos,Fourier,20)(x),20)-Fun(cos,20)) <100eps()
@test norm(Fun(x->Fun(cos,Fourier)(x))-Fun(cos)) <100eps()
@test norm(diff(Fun(cos,Fourier))+Fun(sin,Fourier)) < 100eps()
@test norm(Fun(x->Fun(cos,Laurent)(x))-Fun(cos)) <100eps()
@test norm(diff(Fun(cos,Laurent))+Fun(sin,Laurent)) < 100eps()


@test norm(diff(Fun(cos,Circle()))+Fun(sin,Circle()))<100eps()

f=Fun(exp,Circle());

@test norm(diff(f)-f)<100eps()
@test norm(integrate(f)+1-f)<100eps()

f=Fun(x->exp(-10sin((x-.1)/2)^2),Fourier)
@test_approx_eq real(f)(.1) f(.1)



@test_approx_eq (Fun(z->sin(z)*cos(1/z),Circle())*Fun(z->exp(z)*airyai(1/z),Circle()))(exp(.1im)) (z->sin(z)*cos(1/z)*exp(z)*airyai(1/z))(exp(.1im))

## Calculus

f=Fun(t->cos(t),CosSpace)
D=Derivative(space(f))
@test_approx_eq (D*f)(.1) -sin(.1)
@test_approx_eq diff(f)(.1) -sin(.1)

f=Fun(t->sin(t),SinSpace)
D=Derivative(space(f))
@test_approx_eq (D*f)(.1) cos(.1)
@test_approx_eq diff(f)(.1) cos(.1)

f=Fun(cos,Fourier)
@test norm((Derivative(space(f))^2)*f+f)<10eps()



## Taylor

for d in (Circle(),Circle(0.5),Circle(-0.1,2.))
    S=Taylor(d)
    D=Derivative(S)
    ef=Fun(exp,S)
    @test norm((D*ef-ef).coefficients)<1000eps()
    @test norm((D^2*ef-ef).coefficients)<100000eps()
    u=[Evaluation(S,0.),D-I]\[1.]
    @test norm((u-ef).coefficients)<100eps()
    @test norm((Integral(S)*Fun(exp,S)+ef.coefficients[1]-ef).coefficients)<100eps()


    f=Fun(z->exp(1/z)-1,Hardy{false}(d))
    df=Fun(z->-1/z^2*exp(1/z),Hardy{false}(d))
    @test norm((Derivative()*f-df).coefficients)<1000eps()
    @test norm((Derivative()^2*f-diff(df)).coefficients)<100000eps()
    @test norm((diff(f)-df).coefficients)<1000eps()
end



## Multiplication

s=Fun(t->(sin(t)+sin(2t))*cos(sin(t)),SinSpace)
b=Fun(t->(sin(t)+sin(3t)),SinSpace)
@test_approx_eq (s*s)(.1) s(.1)^2
@test_approx_eq (s*b)(.1) s(.1)*b(.1)

s=Fun(t->(cos(t)+cos(2t))*cos(cos(t)),CosSpace)
b=Fun(t->(1+cos(t)+cos(3t)),CosSpace)

@test_approx_eq (s*s)(.1) s(.1)^2
@test_approx_eq (s*b)(.1) s(.1)*b(.1)

s=Fun(t->(cos(t)+cos(2t))*cos(cos(t)),CosSpace)
b=Fun(t->(sin(t)+sin(3t)),SinSpace)

@test_approx_eq (s*b)(.1) s(.1)*b(.1)


s=Fun(t->(sin(t)+sin(2t))*cos(sin(t)),SinSpace)
b=Fun(t->(1+cos(t)+cos(3t)),CosSpace)

@test_approx_eq (s*b)(.1) s(.1)*b(.1)



a=Fun(t->exp(cos(t)+sin(t)),Fourier)
b=Fun(t->sin(t)+cos(3t)+1,Fourier)

@test_approx_eq (a*b)(.1) a(.1)*b(.1)

a=Fun(t->exp(cos(t)),CosSpace)
b=Fun(t->sin(t)+cos(3t)+1,Fourier)

@test_approx_eq (a*b)(.1) a(.1)*b(.1)

a=Fun(t->sin(sin(t)),SinSpace)
b=Fun(t->sin(t)+cos(3t)+1,Fourier)

@test_approx_eq (a*b)(.1) a(.1)*b(.1)
