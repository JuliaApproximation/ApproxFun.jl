using ApproxFun, Base.Test

@test norm(Fun(x->FFun(cos)[x])-Fun(cos)) <100eps()


@test norm(diff(FFun(cos))+FFun(sin)) < 100eps()

@test norm(diff(FFun(cos,Circle()))+FFun(sin,Circle()))<100eps()

f=Fun(exp,Circle());

@test norm(diff(f)-f)<100eps()
@test norm(integrate(f)+1-f)<100eps()

f=FFun(x->exp(-10sin((x-.1)/2)^2))
@test_approx_eq real(f)[.1] f[.1]



## Calculus

f=Fun(t->cos(t),CosSpace)
D=Derivative(space(f))
@test_approx_eq (D*f)[.1] -sin(.1)
@test_approx_eq diff(f)[.1] -sin(.1)

f=Fun(t->sin(t),SinSpace)
D=Derivative(space(f))
@test_approx_eq (D*f)[.1] cos(.1)
@test_approx_eq diff(f)[.1] cos(.1)

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
end

