using ApproxFun        # solves 2u''' + uu'' = 0 ,  u(0) = u'(0) = 0, u'(∞) = 1
x=Fun(identity,[0.,4π])
d=domain(x)
B=[ldirichlet(d),lneumann(d),rneumann(d)]
D=Derivative(d)
κ = 0.33205733621519630
u0 = (1//2) * κ * x.^2 # us = (1//2) * κ * x.^2 -(1//240)* κ^2 * x.^5 + (11//161280) * κ^3 * x.^8
u = 0.5x^2

f = (u)->(2.0*D^3*u + u*D^2*u)
df = (u)->(2.0*D^3 + u*D^2 + D^2*u)

ApproxFun.plot([u,u0])

for k=1:10
  u  -= [B; df(u)]\[u[0.];diff(u,1)[0.];diff(u,1)[d.b]-1.;f(u)]
end
norm(f(u)) # should be zero
abs(diff(u,2)[0.]-κ) # should be also zero
ApproxFun.plot([u,u0])

m = 0.10                # m ∈ [-0.0905, 2]
β = 2m/(1+m)
F = (u)->(2.0*D^3*u + u*D^2*u - β*(1.0 - (D*u)*(D*u)))
dF = (u)->(2.0*D^3 + u*D^2 + D^2*u + 2β*D*u*D)
v=u
norm(F(v)) # should be non-zero

for k=1:10
  v  -= [B; dF(v)]\[v[0.];diff(v,1);diff(v,1)[d.b]-(1.-exp(-d.b));F(v)]
end
norm(F(v)) # should be zero

ApproxFun.plot([u,v]) # Blasius and Falkner-Skan solutions

# TODO: add MHD variant with additional nonlinearity
