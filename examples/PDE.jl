# # Solving PDEs

# We solve
# ```math
# (Δ + 100)u(x,y) = 0
# ```
# on the rectangle `-1..1 × -1..1`,
# subject to the condition that ``u=1`` on the boundary of the domain.

include("PDE1.jl")

# Plot the solution
import Plots
Plots.plot(u; st=:surface, legend=false)
