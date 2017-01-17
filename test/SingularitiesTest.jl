using ApproxFun, Base.Test


x=Fun(identity);
@test_approx_eq sqrt(cos(π/2*x))(.1) sqrt(cos(.1π/2))


x=Fun(identity,-2..2)
u=sqrt(4-x.^2)/(2π)
@test_approx_eq u(.1) sqrt(4-.1^2)/(2π)
@test_approx_eq sum(u) 1

#this call threw an error, which we check
@test length(values(u))==1


f=Fun(x->x.*cot(π*x/2))
x=Fun(identity)
u=Fun(JacobiWeight(1.,1.,Interval()),(f./(1-x.^2)).coefficients)
@test_approx_eq 1./(.1.*cot(π*.1/2)) (1./u)(.1)

@test_approx_eq (x./u)(.1) tan(π*.1/2)


f=Fun(x->exp(-x.^2),Line(0.,0.,-.5,-.5),400)
@test_approx_eq sum(f) sqrt(π)

f=Fun(x->exp(x)/sqrt(1-x.^2),JacobiWeight(-.5,-.5))
@test_approx_eq f(.1) (x->exp(x)/sqrt(1-x.^2))(.1)



S=JacobiWeight(-1.,-1.,Chebyshev(0..1))
D=Derivative(S)

f=Fun(S,Fun(exp,0..1).coefficients)
x=.1
@test_approx_eq f(x) exp(x)*x^(-1).*(1-x)^(-1)/4
@test_approx_eq (D*f)(x) -exp(x)*(1+(x-3)*x)/(4*(x-1)^2*x^2)


S=JacobiWeight(-1.,0.,Chebyshev(0..1))
D=Derivative(S)

f=Fun(S,Fun(exp,0..1).coefficients)
x=.1
@test_approx_eq f(x) exp(x)*x^(-1)/2
@test_approx_eq (D*f)(x) exp(x)*(x-1)/(2x^2)



## ODEs

## f/g bugs

println("    Jacobi singularity tests")

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


println("    Ray tests")

## Ray

@test Inf in Ray()   # this was a bug

f=Fun(x->exp(-x),0..Inf)
@test_approx_eq f'(.1) -f(.1)

x=Fun(identity,Ray())
f=exp(-x)
u=integrate(f)
@test_approx_eq (u(1.)-u(0)-1) -f(1)



x=Fun(identity,Ray())
f=x^(-0.123)*exp(-x)
@test_approx_eq integrate(f)'(1.) f(1.)


@test_approx_eq_eps sum(Fun(sech,0..Inf)) sum(Fun(sech,0..40)) 1000000eps()


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

x=Fun(identity,-1..1)
f=exp(x+1)-1
@test_approx_eq log(f)(0.1) log(f(0.1))


x=Fun(identity,0..1)
f=exp(x)-1
@test_approx_eq log(f)(0.1) log(f(0.1))


## Test divide sing

x=Fun(identity,0..1)
@test_approx_eq Fun(exp(x)/x-1/x,Chebyshev)(0.1) (exp(0.1)-1)/0.1

x=Fun(identity,0..1)
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
@test_approx_eq linesum(w*log(abs(x-2.)))/π (88.5579588360686)



## Dirac Space

a,b=DiracDelta(0.),DiracDelta(1.)
f=Fun(exp)
g=a+0.2b+f
@test_approx_eq pieces(g)[2](0.) 1.
@test_approx_eq g(.1) exp(.1)
@test_approx_eq sum(g) (sum(f)+1.2)


#Checks prevoius bug
δ=DiracDelta()
x=Fun()
w=sqrt(1-x^2)
w+δ


## PointSpace

f=Fun(x->(x-0.1),ApproxFun.PointSpace([0,0.1,1]))
@test roots(f) == [0.1]

a=Fun(exp,space(f))
@test f/a == Fun(x->(x-0.1)*exp(-x),space(f))

g = f + Fun(2..3)
h = a + Fun(2..3)

@test norm((g/h - ((f/a) + Fun(1,2..3))).coefficients) ≤ 10eps()







## multiplicities
x=Fun(identity,-1..1)
@test_approx_eq (1/x^2)(0.1) 100.
@test_approx_eq (1/x^2)(-0.1) 100.

fc=x*(1+x)^2
@test_approx_eq (1/fc)(0.1) 1/fc(0.1)

fc=x*(1-x)^2
@test_approx_eq (1/fc)(0.1) 1/fc(0.1)

## erf(sqrt(x))

x=Fun(0..1)
@test_approx_eq erf(sqrt(x))(0.1) erf(sqrt(0.1))
@test_approx_eq erfc(sqrt(x))(0.1) erfc(sqrt(0.1))


## norm(u-x)

@test_approx_eq norm(Fun(exp,Legendre(0..1))+sqrt(x)) 2.491141949903508



## Test Jacobi special conversions



S1,S2=JacobiWeight(3.,1.,Jacobi(1.,1.)),JacobiWeight(1.,1.,Jacobi(0.,1.))
f=Fun(S1,[1,2,3.])
C=Conversion(S1,S2)
Cf=C*f
@test_approx_eq Cf(0.1) f(0.1)


S1,S2=JacobiWeight(3.,2.,Jacobi(1.,1.)),JacobiWeight(1.,1.,Jacobi(0.,0.))
f=Fun(S1,[1,2,3.])
C=Conversion(S1,S2)
Cf=C*f
@test_approx_eq Cf(0.1) f(0.1)



## roots of log(abs(x-y))

x=Fun(-2..(-1))
@test_approx_eq roots(abs(x+1.2)) [-1.2]

f=abs(x+1.2)

@test norm(abs(f)-f)<10eps()
@test norm(sign(f)-Fun(1,space(f)))<10eps()


@test_approx_eq log(f)(-1.3) log(abs(-1.3+1.2))
@test_approx_eq log(f)(-1.1) log(abs(-1.1+1.2))
