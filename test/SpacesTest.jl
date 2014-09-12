using ApproxFun, Base.Test

import ApproxFun: ChebyshevDirichletSpace,UltrasphericalSpace


@test_approx_eq IFun(exp,ChebyshevDirichletSpace{1,1})[.1] exp(.1)
@test_approx_eq IFun(IFun(exp,ChebyshevDirichletSpace{1,1}),UltrasphericalSpace{1})[.1] exp(.1)

d=Interval()
sp=ChebyshevDirichletSpace{1,1}(d)
B=dirichlet(sp)
D=diff(d)
L=D^2+I

norm(([B,L]\[1.])-([dirichlet(d),L]\[1.]))