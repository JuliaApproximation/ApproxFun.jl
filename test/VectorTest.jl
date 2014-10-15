using ApproxFun,Base.Test

import ApproxFun:interlace


d=Interval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();
A=[B 0;
   0 B;
   D-I 2.I;
   0 D+I];
Ai=interlace(A);
f=Fun(x->[exp(x),cos(x)],d)

u=Ai\{0.,0.,f};
u1=vec(u)[1];u2=vec(u)[2];
f1=vec(f)[1];f2=vec(f)[2];
@test norm(diff(u2)+u2-f2)<eps()
@test norm(diff(u1)-u1+2.u2-f1)<eps()




## Matrix exponential

n=4
d=fill(Interval(0.,1.),n)
B=Evaluation(d,0.)
D=Derivative(d)
A=rand(n,n)
L=[B;D-A]
u=L\eye(n)
@test norm(evaluate(u,1.)-expm(A))<eps(1000.)



