using Plots,ApproxFun
# The Blasius boundary layer problem (https://en.wikipedia.org/wiki/Blasius_boundary_layer)
# describes the boundary layer on a semi-infinite plane with a head-on steady fluid flow.
# When properly non-dimensionalized, it presents itself as a relatively simple nonlinear BVP.
# The Falkner-Skan problem generalizes the Blasius boundary layer by considering a wedge at an
# angle to the steady flow.
# In both cases, no-slip conditions are in place at the wall and the boundary layer
# grows linearly (u'(∞)→1) in the far-field.
# There are more variants which include different nonlinearities which arise in magneto-
# hydrodynamics, for example.

# We begin with Blasius:
#   solves 2u''' + uu'' = 0 ,  u(0.) = u'(0.) = 0, u'[∞] = 1
x=Fun(identity,0..4π)
d=domain(x)
B=[ldirichlet(d),lneumann(d),rneumann(d)]
D=Derivative(d)
κ = 0.33205733621519630     # This is diff(u,2)(0.), due to Boyd.
u0 = (1//2) * κ * x^2      # First-order approximation
u = 0.5x^2                  # Other initial solutions may fail to converge

f = (u)->(2.0*D^3*u + u*D^2*u)
df = (u)->(2.0*D^3 + u*D^2 + D^2*u)

plot([u,u0])      # Initial estimate and first-order approximation

for k=1:10
  u  -= [B; df(u)]\[u(0.);u'(0.);u'(d.b)-1.;f(u)]
end
norm(f(u))             # should be zero
abs(u''(0.)-κ)         # should also be zero
plot([u,u0])           # Solution and First-order approximation
# Now for Falkner-Skan:
#     solves 2u''' + uu'' - β(1 - (u')^2) = 0 ,  u(0.) = u'(0.) = 0, u'[∞] = 1
m = 0.11                # m ∈ [-0.0905, 2]
β = 2m/(1+m)
F = (u)->(2.0*D^3*u + u*D^2*u - β*(1.0 - (D*u)*(D*u)))
dF = (u)->(2.0*D^3 + u*D^2 + D^2*u + 2β*D*u*D)
v=u
norm(F(v))               # should be non-zero, as we've perturbed the Blasius equation

for k=1:10
  v  -= [B; dF(v)]\[v(0.);v'(0.);v'(d.b)-1.;F(v)]
end
norm(F(v))               # should be zero

plot([u,v])    # Blasius and Falkner-Skan solutions
