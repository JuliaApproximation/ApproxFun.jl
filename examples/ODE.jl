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

include("ODE_BVP.jl")

# We plot the solution
import Plots
Plots.plot(u, xlabel="x", ylabel="u(x)", legend=false)


# ## Increasing Precision

# Solving differential equations with high precision types is available.
# The following calculates ``e`` to 300 digits by solving the ODE ``u^\prime = u``:

include("ODE_increaseprec.jl")

Plots.plot(u; legend=false, xlabel="x", ylabel="u(x)")
