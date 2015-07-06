using ApproxFun, Base.Test

import ApproxFun: ChebyshevDirichlet,Ultraspherical,space


@test_approx_eq Fun(exp,ChebyshevDirichlet{1,1})[.1] exp(.1)
@test_approx_eq Fun(Fun(exp,ChebyshevDirichlet{1,1}),Ultraspherical{1})[.1] exp(.1)

d=Interval()
sp=ChebyshevDirichlet{1,1}(d)
B=dirichlet(sp)
D=diff(d)
L=D^2+I

@test norm(([B;L]\[1.])-([dirichlet(d);L]\[1.])) <eps()

f=Fun(t->cos(t)+cos(3t),CosSpace)

@test (f.*f-Fun(t->(cos(t)+cos(3t))^2,CosSpace)).coefficients|>norm <100eps()



f=Fun(exp,Taylor(Circle()))
g=Fun(z->1./(z-.1),Hardy{false}(Circle()))
@test_approx_eq (f[1.]+g[1.]) (exp(1.) + 1./(1-.1))


## Periodic
f=Fun(x->exp(-10sin((x-.1)/2)^2),Laurent)
@test_approx_eq f[.5] (Conversion(space(f),Fourier(domain(f)))*f)[.5]
@test_approx_eq f[.5] Fun(f,Fourier)[.5]



## PiecewiseSPace

x=Fun(identity,[-10.,0.,1.,15.])
sp=space(x)
D=Derivative(sp)

u=[dirichlet(sp);
    D^2-x]\[airyai(-10.)];

@test_approx_eq u[0.] airyai(0.)

s=Fun(sin,[-2.,2.])|>abs
c=Fun(cos,[-2.,2.])|>abs
sc=Fun(x->abs(sin(x))+abs(cos(x)),[-2,-π/2,0,π/2,2])
@test norm(sc-(c+s))<100eps()




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
D=Derivative(S)
u=[ldirichlet(S),D-I]\[exp(1.)]


@test_approx_eq u[1.1] exp(1.1)
@test_approx_eq u[3.4] exp(3.4)
@test_approx_eq last(u) exp(4)


d=PiecewiseInterval(0,1.,1.+im,im,0.)
@test_approx_eq Fun(exp,d)[.1] exp(.1)





## Triple SumSpace

x=Fun()
w=log(1-x)+sqrt(1-x^2)
f=w+x
@test_approx_eq f[0.1] (w[0.1]+0.1)
@test_approx_eq (w+1)[0.1] (w[0.1]+1)
@test_approx_eq (w+x+1)[0.1] (w[0.1]+1.1)
@test_approx_eq ((w+x)+1)[0.1] (w[0.1]+1.1)






## Piecewise + Cosntant

Γ=Circle()∪Circle(0.0,0.4)
G=Fun(z->in(z,Γ[2])?[1 0; -1/z 1]:[z 0; 0 1/z],Γ)   # Before the 80 wasn’t specified causing inconsistency
@test_approx_eq (G-I)[1.] (G[1.]-I)
