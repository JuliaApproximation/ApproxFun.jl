# # Nonlinear Boundary Value Problem

# There is preliminary support for nonlinear equations, via Newton iteration in function space.
# Here is a simple two-point boundary value problem:

# We solve
# ```math
# Du = 0.001u^{\prime\prime} + 6(1-x^2)u^{\prime} + u^2 - 1 = 0,
# ```
# subject to the boundary conditions ``u(-1)=1`` and ``u(1)=-0.5``.

include("NonlinearBVP1.jl")

# We plot the solution
import Plots
Plots.plot(u; title = "Solution", xlabel="x", ylabel="u(x)", legend=false)

# ## System of nonlinear differential equations
# One can also solve a system of nonlinear ODEs with potentially nonlinear boundary conditions:
N(u1, u2) = [u1'(0) - 0.5*u1(0)*u2(0);
                u2'(0) + 1;
                u1(1) - 1;
                u2(1) - 1;
                u1'' + u1*u2;
                u2'' - u1*u2]

function nbvpsolver2()
    x = Fun(0..1)
    u10 = one(x)
    u20 = one(x)
    newton(N, [u10,u20])
end
u1,u2 = nbvpsolver2();

# Plot the solutions
Plots.plot(u1, label="u1", xlabel="x")
Plots.plot!(u2, label="u2")
