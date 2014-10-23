using ApproxFun

x=Fun(x->x)
d=domain(x)
D=diff(d)

u=[dirichlet(d),
   1/70*D^2-x*D+I] \ [1.,2.]

plot(u)