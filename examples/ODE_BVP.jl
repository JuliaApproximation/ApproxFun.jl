using ApproxFun
using SpecialFunctions
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
B_vals = [airyai(a), airyai(b)];

# The main step, where we solve the differential equation
u = [B; L] \ [B_vals, 0];

#src we compare the result with the expected solution (which is an Airy function)
using Test #src
@test ≈(u(0), airyai(0); atol=10000eps()) #src
