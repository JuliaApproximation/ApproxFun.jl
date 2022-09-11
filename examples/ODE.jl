# # Solving ODEs

# ## Boundary value problem
# We solve the [Airy equation](https://en.wikipedia.org/wiki/Airy_function)
# ```math
# \frac{d^2 y}{dx^2}-x y = 0,
# ```
# in an interval `a..b`, subject to the boundary conditions ``y(a)=\mathrm{Ai}(a)``
# and ``y(b)=\mathrm{Ai}(b)``, where ``\mathrm{Ai}`` represents the Airy function of the
# first kind.

using ApproxFun

# Construct the domain
a, b = -20, 10
d = a..b;

# We construct the Airy differential operator ``L = d^2/dx^2 - x``:
x = Fun(d);
D = Derivative(d);
L = D^2 - x;

# We impose boundary conditions by setting the function at
# the boundaries equal to the values of the Airy function at these points.
# First, we construct the Dirichlet boundary condition operator,
# that evaluates the function at the boundaries
B = Dirichlet(d);

# The right hand side for the boundary condition is obtained by evaluating the Airy function
# We use `airyai` from `SpecialFunctions.jl` for this
using SpecialFunctions
B_vals = [airyai(a), airyai(b)];

# The main step, where we solve the differential equation
u = [B; L] \ [B_vals, 0];

# We plot the solution
import Plots
Plots.plot(u, xlabel="x", ylabel="u(x)", legend=false)

#src # Finally, we compare the result with the expected solution (which is an Airy function)
using Test #src
@test ≈(u(0), airyai(0); atol=10000eps()) #src

# ## Increasing Precision

# Solving differential equations with high precision types is available.
# The following calculates ``e`` to 300 digits by solving the ODE ``u^\prime = u``:

u = setprecision(1000) do
    d = BigFloat(0)..BigFloat(1)
    D = Derivative(d)
    [ldirichlet(); D-1] \ [1; 0]
end
Plots.plot(u; legend=false, xlabel="x", ylabel="u(x)")

@test u(1) ≈ exp(BigFloat(1)) #src
