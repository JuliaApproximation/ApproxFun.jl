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

include("NonlinearBVP2.jl")

# Plot the solutions
import Plots
Plots.plot(u1, label="u1", xlabel="x")
Plots.plot!(u2, label="u2")
