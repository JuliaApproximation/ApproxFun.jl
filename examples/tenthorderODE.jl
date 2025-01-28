using Plots,ApproxFun

x=Fun()
D=Derivative()


L=D^10 + cosh(x)*D^8 + x^3*D^6 + x^4*D^4 + cos(x)*D^2 + x^2
d=Interval()
B=[dirichlet(d);neumann(d);[Evaluation(Interval(),false,k) for k=2:4];
           [Evaluation(Interval(),true,k) for k=2:4]]

u=[B;L]\[0.;0.;1.;1.;zeros(6);exp(x)]


ApproxFun.plot(u)
