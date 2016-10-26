using ApproxFun, Base.Test
    import ApproxFun: ChebyshevDirichlet,Ultraspherical,space,testspace,testbandedoperator,testcalculus,testtransforms


testtransforms(ChebyshevDirichlet{1,1}())

@test_approx_eq Fun(exp,ChebyshevDirichlet{1,1})(.1) exp(.1)
@test_approx_eq Fun(Fun(exp,ChebyshevDirichlet{1,1}),Ultraspherical(1))(.1) exp(.1)

d=Interval()
sp=ChebyshevDirichlet{1,1}(d)
B=dirichlet(sp)
D=Derivative(d)
L=D^2+I

@test_approx_eq B[1][1:3] [1.;-1.;0.]
@test_approx_eq B[2][1:3] [1.;1.;0.]
@test_approx_eq B[2][1:1,1:3] [1. 1. 0.]


@test_approx_eq csc(2)sin(1 - 0.1)  ([dirichlet(d);L]\[1.])(0.1)
@test_approx_eq csc(2)sin(1 - 0.1)  ([B;L]\[1.])(0.1)

@test norm(([B;L]\[1.])-([dirichlet(d);L]\[1.])) <10eps()




## PiecewiseSPace

x=Fun(identity,[-1.,0.,1.])
sp=space(x)
testtransforms(sp;minpoints=2)

D=Derivative(sp)
testbandedoperator(D)

u=[dirichlet(sp);
    D^2]\[1];
u2=[dirichlet();Derivative(Chebyshev())^2]\[1.]
@test_approx_eq u(0.) u2(0.)

x=Fun(identity,[-10.,0.,1.,15.])
sp=space(x)
D=Derivative(sp)

u=[dirichlet(sp);
    D^2-x]\[airyai(-10.)];

@test_approx_eq u(0.) airyai(0.)

s=Fun(sin,[-2.,2.])|>abs
c=Fun(cos,[-2.,2.])|>abs
sc=Fun(x->abs(sin(x))+abs(cos(x)),[-2,-π/2,0,π/2,2])
@test norm(sc-(c+s))<100eps()

# max/min creates breakpoints

x=Fun()
g=4*(x-.2)
f=max(-1,g)
f2=min(f,1)
f3=Fun(x->x<-0.05?-1.0:(x<0.45?4*(x-.2):1),[-1.0;-0.05;0.45;1.0])
@test norm(f2(collect(linspace(-1,1,10)))-f3(collect(linspace(-1,1,10)))) < 2eps()

x=Fun(identity,[im,0.,1.])
@test_approx_eq x(0.5) 0.5
@test_approx_eq x(0.5im) 0.5im

@test Fun(Fun(1.0),space(x))(0.5) == 1.0
@test Fun(Fun(1.0),space(x))(0.5im) == 1.0

@test_approx_eq (x+1)(0.5) 1.5
@test_approx_eq (x-1)(0.5) -0.5
@test_approx_eq (1-x)(0.5) 0.5



@test_approx_eq sqrt(1-x)(0.2im) sqrt(1-0.2im)
@test_approx_eq sqrt(1-x)(0.2) sqrt(1-0.2)

w=2/(sqrt(1-x)*sqrt(1+im*x))
@test_approx_eq w(0.2im) 2/(sqrt(1-0.2im)*sqrt(1+im*(0.2im)))
@test_approx_eq w(0.2) 2/(sqrt(1-0.2)*sqrt(1+im*(0.2)))

## Equivalent spaces

@test norm(Fun(cos,Chebyshev)-Fun(cos,Jacobi(-0.5,-0.5)))<100eps()
@test norm(Fun(cos,Chebyshev)-Fun(cos,JacobiWeight(0,0)))<100eps()
@test norm(Fun(cos,Jacobi(-0.5,-0.5))-Fun(cos,JacobiWeight(0,0))) < 100eps()
@test norm(Fun(cos,Chebyshev)-Fun(cos,JacobiWeight(0,0,Jacobi(-0.5,-0.5))))<100eps()
@test norm(Fun(cos,Jacobi(-0.5,-0.5))-Fun(cos,JacobiWeight(0,0,Jacobi(-0.5,-0.5))))<100eps()



## ContinuousSpace

import ApproxFun: PiecewiseInterval,ContinuousSpace

d=PiecewiseInterval(1.,2.,3.,4.)
S=ContinuousSpace(d)
testtransforms(S;minpoints=3,invertibletransform=false)

D=Derivative(S)
testbandedoperator(D)

u=[ldirichlet(S),D-I]\[exp(1.)]


@test_approx_eq u(1.1) exp(1.1)
@test_approx_eq u(3.4) exp(3.4)
@test_approx_eq last(u) exp(4)


d=PiecewiseInterval(0,1.,1.+im,im,0.)
@test_approx_eq Fun(exp,d)(.1) exp(.1)





## Triple SumSpace

x=Fun()
w=log(1-x)+sqrt(1-x^2)
f=w+x
@test_approx_eq f(0.1) (w(0.1)+0.1)
@test_approx_eq (w+1)(0.1) (w(0.1)+1)
@test_approx_eq (w+x+1)(0.1) (w(0.1)+1.1)
@test_approx_eq ((w+x)+1)(0.1) (w(0.1)+1.1)


## SumSpace bug

dsp=JacobiWeight(1.,0.,Jacobi(0.,1.,[0.,1.]))⊕JacobiWeight(0.5,0.,Jacobi(-0.5,0.5,[0.,1.]))
rsp=Legendre([0.,1.])⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,[0.,1.]))

testcalculus(sp)

C=Conversion(dsp,rsp)
f=Fun([1.,2.,3.,4.,5.],dsp)
@test_approx_eq f(0.1) (C*f)(0.1)






## Piecewise + Cosntant
using Base.Test

Γ=Circle()∪Circle(0.0,0.4)
o=ones(Γ)
@test_approx_eq o(1.) 1.0
@test_approx_eq o(0.4) 1.0

G=Fun(z->in(z,Γ[2])?[1 0; -1/z 1]:[z 0; 0 1/z],Γ)
@test_approx_eq (G-I)(1.) (G(1.)-I)


## Previoius seffdault

x=Fun(identity,[-1.,1.])
f=x+sin(2x)*sqrt(1-x^2)
@test_approx_eq f(0.1) 0.1+sin(2*0.1)*sqrt(1-0.1^2)


## Check multiple piecewisesapce

x=Fun(identity,[-3,-2])+Fun(identity,[2,3])
w=sqrt(9-x^2)
f=w+Fun()
@test_approx_eq (f+w)(2.5) 2w(2.5)
@test_approx_eq (f+w)(.5) f(.5)



## Check Jacobi recurrence bug

S=Jacobi(-.5,.5)
f=Fun(exp,S)
@test_approx_eq f(0.1) exp(0.1)


## Check cancel conversion works
x=Fun([0.,1.])
f=exp(x)-1
Fun(f,JacobiWeight(1.,0.,[0.,1.]))


## Hermite
f=Fun(x->x+x^2,Hermite())
@test_approx_eq f(1.) 2.



## Arc exp

z=Fun(identity,Arc(0.,.1,0.,π/2))
@test_approx_eq exp(z)(0.1exp(0.2im)) exp(0.1exp(0.2im))



## Extending function

Γ=Interval(-im,1.0-im)∪Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6)))∪Circle(2.0,0.2)

@test isempty(Γ[1]\Γ[1])
@test Γ\Γ[1] == Γ[2]∪Γ[3]

@test norm(Fun(ones(Γ[1]),Γ) - Fun(x->x ∈ Γ[1]?1.0:0.0,Γ)) == 0


## Line

f=Fun(z->2exp(z^2),PeriodicLine(0.,π/2))
@test_approx_eq f(1.1im) 2exp(-1.1^2)


f=Fun(z->2exp(z^2),Line(0.,π/2))
@test_approx_eq f(1.1im) 2exp(-1.1^2)



## Exp for Γ

a=1+10*im;b=2-6*im
d=Curve(Fun(x->1+a*x+x^2+b*x^3))

x=Fun(d)

@test_approx_eq exp(x)(1+a*0.1+0.1^2+b*0.1^3) exp(1+a*0.1+0.1^2+b*0.1^3)


## ChebyshevDirichlet multiplication

S=ChebyshevDirichlet()
x=Fun()
@test norm((ApproxFun.Recurrence(S)*Fun(exp,S)-Fun(x->x*exp(x),S)).coefficients) < 100eps()
@test norm((x*Fun(exp,S)-Fun(x->x*exp(x),S)).coefficients) < 100eps()
