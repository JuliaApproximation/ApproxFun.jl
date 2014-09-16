using ApproxFun

# solves u" + 2(1-x^2)u + u^2 = 1 ,  u(-1) = u(1) = 0

x=Fun(identity)
d=domain(x)
B=dirichlet(d)
D=diff(d)

# Solves Lu + g(u)-1==0

L=D^2 + 2(1-x.^2)*D
g=u->u.^2;gp=u->2u

u=0.x   #initial guess is zero

for k=1:5
        u=u-[B,L+gp(u)]\[0.,0.,L*u+g(u)-1.];
end

norm(diff(u,2) + 2(1-x.^2).*diff(u) + g(u) -1)  # This equals 0.0