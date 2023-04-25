# # Solving ODEs

# ## Boundary value problem
# We solve the [Airy equation](https://en.wikipedia.org/wiki/Airy_function)
# ```math
# \frac{d^2 y}{dx^2}-x y = 0,
# ```
# in an interval `a..b`, subject to the boundary conditions ``y(a)=\mathrm{Ai}(a)``
# and ``y(b)=\mathrm{Ai}(b)``, where ``\mathrm{Ai}`` represents the Airy function of the
# first kind.

include("ODE_BVP.jl")

# We plot the solution
import Plots
Plots.plot(u, xlabel="x", ylabel="u(x)", legend=false)


# ## Initial value problem
# ### Undamped harmonic oscillator
# ```math
# \frac{d^2 u}{dt^2} + 4 u = 0,
# ```
# in an interval `0..2pi`, subject to the initial conditions ``u(0)=0``
# and ``u'(0)=2``.

include("ODE_IVP.jl")

# We plot the solution
t = a:0.1:b
Plots.plot(t, u.(t), xlabel="t", label="u(t)")
Plots.plot!(t, sin.(2t), label="Analytical", seriestype=:scatter)

# ### Damped harmonic oscillator
# ```math
# \frac{d^2 y}{dt^2} + 2\zeta\omega_0\frac{dy}{dt} + \omega_0^2 y = 0,
# ```
# from ``t=0`` to ``t=T``, subject to the initial conditions ``y(0)=4``
# and ``y'(0)=0``.

include("ODE_IVP2.jl")

# Plot the solution
import Plots
Plots.plot(y, xlabel="t", ylabel="y(t)", legend=false)

# ## Increasing Precision

# Solving differential equations with high precision types is available.
# The following calculates ``e`` to 300 digits by solving the ODE ``u^\prime = u``:

include("ODE_increaseprec.jl")
import Plots
Plots.plot(u; legend=false, xlabel="x", ylabel="u(x)")
