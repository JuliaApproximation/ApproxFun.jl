x=Fun(x->x);
d=x.domain;
e=1./70;
I=eye(d);
D=diff(d);
D2=diff(d,2);
A=e*D2-x*D+I;
B=dirichlet(d);

u=[B,A] \ [1.,2.,0.];

plot(u)