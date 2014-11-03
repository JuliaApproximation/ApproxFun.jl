using ApproxFun

x=Fun(identity)
d=domain(x)
D=Derivative(d)

u=[dirichlet(d),
   1/70*D^2-x*D+I] \ [1.,2.]

ApproxFun.plot(u)