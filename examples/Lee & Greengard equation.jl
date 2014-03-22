using ApproxFun

x=Fun(x->x);
d=x.domain;
D=diff(d);

u=[dirichlet(d),
  1./70*D^2-x*D+1.] \ [1.,2.];

plot(u)