using ApproxFun, Base.Test


x=Fun(identity);
@test_approx_eq sqrt(cos(π/2*x))[.1] sqrt(cos(.1π/2))


x=Fun(identity,[-2.,2.])
u=sqrt(4-x.^2)/(2π)
@test_approx_eq u[.1] sqrt(4-.1^2)/(2π)
@test_approx_eq sum(u) 1

values(u)


f=Fun(x->x.*cot(π*x/2))
x=Fun(identity)
u=Fun((f./(1-x.^2)).coefficients,JacobiWeight(1.,1.,Interval()))
@test_approx_eq 1./(.1.*cot(π*.1/2)) (1./u)[.1]

@test_approx_eq (x./u)[.1] tan(π*.1/2)


f=Fun(x->exp(-x.^2),Line(0.,0.,-.5,-.5),400)
@test_approx_eq sum(f) sqrt(π)

f=Fun(x->exp(x)/sqrt(1-x.^2),JacobiWeight(-.5,-.5))
@test_approx_eq f[.1] (x->exp(x)/sqrt(1-x.^2))(.1)



S=JacobiWeight(-1.,-1.,Chebyshev([0.,1.]))
D=Derivative(S)

f=Fun(Fun(exp,[0.,1.]).coefficients,S)
x=.1
@test_approx_eq f[x] exp(x)*x^(-1).*(1-x)^(-1)/4
@test_approx_eq (D*f)[x] -exp(x)*(1+(x-3)*x)/(4*(x-1)^2*x^2)


S=JacobiWeight(-1.,0.,Chebyshev([0.,1.]))
D=Derivative(S)

f=Fun(Fun(exp,[0.,1.]).coefficients,S)
x=.1
@test_approx_eq f[x] exp(x)*x^(-1)/2
@test_approx_eq (D*f)[x] exp(x)*(x-1)/(2x^2)



## ODEs


for ν in (1.,.123,2.,3.5)
    S=JacobiWeight(-ν,0.,Chebyshev([0.,1.]))
    D=Derivative(S)
    x=Fun(identity,domain(S))
    L=(x^2)*D^2+x*D+(x^2-ν^2);
    u=[rdirichlet(S);rneumann(S);L]\[bessely(ν,1.),.5*(bessely(ν-1.,1.)-bessely(ν+1.,1.))]
    @test_approx_eq_eps u[.1] bessely(ν,.1) eps(10000.)*max(abs(u[.1]),1)
    u=[rdirichlet(S),rneumann(S),L]\[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.))]
    @test_approx_eq_eps u[.1] besselj(ν,.1) eps(10000.)*max(abs(u[.1]),1)
end






## f/g bugs

x = Fun(identity)
f = exp(x)./(1-x.^2)

@test_approx_eq f[.1] exp(.1)./(1-.1^2)
f = exp(x)./(1-x.^2).^1
@test_approx_eq f[.1] exp(.1)./(1-.1^2)
f = exp(x)./(1-x.^2).^1.0
@test_approx_eq f[.1] exp(.1)./(1-.1^2)



## 1/f with poles

x=Fun(identity)
f=sin(10x)
g=1/f

@test_approx_eq g[.123] csc(10*.123)





## Ray

f=Fun(x->exp(-x),[0,Inf])
@test_approx_eq diff(f)[.1] -f[.1]

x=Fun(identity,Ray())
f=exp(-x)
u=integrate(f)
@test_approx_eq (u[1.]-u[0]-1) -f[1]



x=Fun(identity,Ray())
f=x^(-0.123)*exp(-x)
@test_approx_eq diff(integrate(f))[1.] f[1.]


@test_approx_eq_eps sum(Fun(sech,[0,Inf])) sum(Fun(sech,[0,40.])) 100000eps()



## Line

f=Fun(x->exp(-x^2),Line())

@test_approx_eq f'[0.1] -2*0.1exp(-0.1^2)
@test_approx_eq (Derivative()*f)[0.1] -2*0.1exp(-0.1^2)
