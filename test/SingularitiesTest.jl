using ApproxFun, Base.Test


x=Fun(identity);
@test_approx_eq sqrt(cos(π/2*x))(.1) sqrt(cos(.1π/2))


x=Fun(identity,[-2.,2.])
u=sqrt(4-x.^2)/(2π)
@test_approx_eq u(.1) sqrt(4-.1^2)/(2π)
@test_approx_eq sum(u) 1

values(u)


f=Fun(x->x.*cot(π*x/2))
x=Fun(identity)
u=Fun((f./(1-x.^2)).coefficients,JacobiWeight(1.,1.,Interval()))
@test_approx_eq 1./(.1.*cot(π*.1/2)) (1./u)(.1)

@test_approx_eq (x./u)(.1) tan(π*.1/2)


f=Fun(x->exp(-x.^2),Line(0.,0.,-.5,-.5),400)
@test_approx_eq sum(f) sqrt(π)

f=Fun(x->exp(x)/sqrt(1-x.^2),JacobiWeight(-.5,-.5))
@test_approx_eq f(.1) (x->exp(x)/sqrt(1-x.^2))(.1)



S=JacobiWeight(-1.,-1.,Chebyshev([0.,1.]))
D=Derivative(S)

f=Fun(Fun(exp,[0.,1.]).coefficients,S)
x=.1
@test_approx_eq f(x) exp(x)*x^(-1).*(1-x)^(-1)/4
@test_approx_eq (D*f)(x) -exp(x)*(1+(x-3)*x)/(4*(x-1)^2*x^2)


S=JacobiWeight(-1.,0.,Chebyshev([0.,1.]))
D=Derivative(S)

f=Fun(Fun(exp,[0.,1.]).coefficients,S)
x=.1
@test_approx_eq f(x) exp(x)*x^(-1)/2
@test_approx_eq (D*f)(x) exp(x)*(x-1)/(2x^2)



## ODEs


for ν in (1.,.123,2.,3.5)
    S=JacobiWeight(-ν,0.,Chebyshev([0.,1.]))
    D=Derivative(S)
    x=Fun(identity,domain(S))
    L=(x^2)*D^2+x*D+(x^2-ν^2);
    u=[rdirichlet(S);rneumann(S);L]\[bessely(ν,1.),.5*(bessely(ν-1.,1.)-bessely(ν+1.,1.))]
    @test_approx_eq_eps u(.1) bessely(ν,.1) eps(10000.)*max(abs(u(.1)),1)
    u=[rdirichlet(S),rneumann(S),L]\[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.))]
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(10000.)*max(abs(u(.1)),1)
end






## f/g bugs

x = Fun(identity)
f = exp(x)./(1-x.^2)

@test_approx_eq f(.1) exp(.1)./(1-.1^2)
f = exp(x)./(1-x.^2).^1
@test_approx_eq f(.1) exp(.1)./(1-.1^2)
f = exp(x)./(1-x.^2).^1.0
@test_approx_eq f(.1) exp(.1)./(1-.1^2)



## 1/f with poles

x=Fun(identity)
f=sin(10x)
g=1/f

@test_approx_eq g(.123) csc(10*.123)





## Ray

f=Fun(x->exp(-x),[0,Inf])
@test_approx_eq diff(f)(.1) -f(.1)

x=Fun(identity,Ray())
f=exp(-x)
u=integrate(f)
@test_approx_eq (u(1.)-u(0)-1) -f(1)



x=Fun(identity,Ray())
f=x^(-0.123)*exp(-x)
@test_approx_eq diff(integrate(f))(1.) f(1.)


@test_approx_eq_eps sum(Fun(sech,[0,Inf])) sum(Fun(sech,[0,40.])) 1000000eps()


#Ei (Exp Integral)

y=Fun(Ray())
q=integrate(exp(-y)/y)
@test_approx_eq (q-last(q))(2.) (-0.04890051070806113)



## Line

f=Fun(x->exp(-x^2),Line())

@test_approx_eq f'(0.1) -2*0.1exp(-0.1^2)
@test_approx_eq (Derivative()*f)(0.1) -2*0.1exp(-0.1^2)




## PeriodicLine

d=PeriodicLine()
D=Derivative(d)
f=Fun(x->sech(x-.1),d)


@test_approx_eq_eps (D*f)(.2) -0.0991717226583897  100000eps()
@test_approx_eq_eps (D^2*f)(.2) -0.9752522555114987  1000000eps()



## LogWeight

x=Fun(identity,[-1.,1.])
f=exp(x+1)-1
@test_approx_eq log(f)(0.1) log(f(0.1))


x=Fun(identity,[0.,1.])
f=exp(x)-1
@test_approx_eq log(f)(0.1) log(f(0.1))


## Test divide sing

x=Fun(identity,[0,1])
@test_approx_eq Fun(exp(x)/x-1/x,Chebyshev)(0.1) (exp(0.1)-1)/0.1

x=Fun(identity,[0,1])
f=1/x
p=integrate(f)
@test_approx_eq (p-p(1.))(0.5) log(0.5)

f=1/(1-x)
p=integrate(f)
@test_approx_eq (p-p(0.))(0.5) -log(1-0.5)



y=Fun(Ray())
@test_approx_eq (y^2)(10.) y(10.)^2
@test_approx_eq 1/y(10.) (1/y)(10.)
@test_approx_eq (1/y^2)(10.) 1/y(10.)^2
@test_approx_eq (-1/y^2)'(10.) 2/(y(10.)^3)
@test_approx_eq exp(-1/y^2)(5.) exp(-1/y(5.)^2)



# catch bug from SIE

a=1+10*im;b=2-6*im
d=Curve(Fun(x->1+a*x+b*x^2))


x=Fun(d)
w=sqrt(abs(first(d)-x))*sqrt(abs(last(d)-x))

@test_approx_eq sum(w/(x-2.))/(2π*im) (-4.722196879007759+2.347910413861846im)