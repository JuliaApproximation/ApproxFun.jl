using ApproxFun

x=Fun(identity)
D=Derivative()

u=[dirichlet();
   1/70*D^2-x*D+I] \ [1,2,0]
