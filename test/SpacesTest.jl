using ApproxFun, Base.Test

import ApproxFun: ChebyshevDirichletSpace,Ultraspherical,space


@test_approx_eq Fun(exp,ChebyshevDirichletSpace{1,1})[.1] exp(.1)
@test_approx_eq Fun(Fun(exp,ChebyshevDirichletSpace{1,1}),Ultraspherical{1})[.1] exp(.1)

d=Interval()
sp=ChebyshevDirichletSpace{1,1}(d)
B=dirichlet(sp)
D=diff(d)
L=D^2+I

@test norm(([B,L]\[1.])-([dirichlet(d),L]\[1.])) <eps()

f=Fun(t->cos(t)+cos(3t),CosSpace)

@test (f.*f-Fun(t->(cos(t)+cos(3t))^2,CosSpace)).coefficients|>norm <100eps()



f=Fun(exp,TaylorSpace(Circle()))
g=Fun(z->1./(z-.1),HardySpace{false}(Circle()))
@test_approx_eq (f[1.]+g[1.]) (exp(1.) + 1./(1-.1))


## Periodic
f=FFun(x->exp(-10sin((x-.1)/2)^2))
@test_approx_eq f[.5] (Conversion(space(f),FourierSpace(domain(f)))*f)[.5]
@test_approx_eq f[.5] Fun(f,FourierSpace)[.5]



## PiecewiseSPace

x=Fun(identity,[-10.,0.,1.,15.])
sp=space(x)
D=Derivative(sp)

u=[dirichlet(sp),
    D^2-x]\[airyai(-10.)];
    
@test_approx_eq u[0.] airyai(0.)

s=Fun(sin,[-2.,2.])|>abs
c=Fun(cos,[-2.,2.])|>abs
@test norm(Fun(x->abs(sin(x))+abs(cos(x)),[-2,-π/2,0,π/2,2])-(c+s))<100eps()

