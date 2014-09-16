using ApproxFun, Base.Test

import ApproxFun: ChebyshevDirichletSpace,UltrasphericalSpace


@test_approx_eq Fun(exp,ChebyshevDirichletSpace{1,1})[.1] exp(.1)
@test_approx_eq Fun(Fun(exp,ChebyshevDirichletSpace{1,1}),UltrasphericalSpace{1})[.1] exp(.1)

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