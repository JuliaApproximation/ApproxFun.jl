# # Nonlinear Boundary Value Problem

# ## Non-linear differential equation
# We solve
# ```math
# Du = 0.001u^{\prime\prime} + 6(1-x^2)u^{\prime} + u^2 - 1 = 0,
# ```
# subject to the boundary conditions ``u(-1)=1`` and ``u(1)=-0.5``.
# We collectively express the system as
# ```math
    # Nu = [u(-1)-1,\,u(1)+0.5,\,Du] = [0,0,0].
# ```

using ApproxFun
using LinearAlgebra

# Define the vector that collates the differential equation and
# the boundary conditions
N(u, x = Fun()) = [u(-1.)-1., u(1.)+0.5, 0.001u'' + 6(1-x^2)u' + u^2 - 1];

# Solve the equation using Newton iteration
function nbvpsolver()
    x = Fun()
    u0 = 0 * x # starting value

    newton(N, u0)
end

u = nbvpsolver();

# We plot the solution
import Plots
Plots.plot(u; title = "Solution", xlabel="x", ylabel="u(x)", legend=false)

#src # We verify that the solution satisfies the differential equation and the boundary conditions
using Test #src
@test norm(N(u)) ≤ 1000eps() #src

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
