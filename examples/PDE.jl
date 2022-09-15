# # Solving PDEs

include("PDE_Poisson.jl")

import Plots
Plots.plot(u; st=:surface, legend=false)

# Many PDEs have weak singularities at the corners,
# in which case it is beneficial to specify a `tolerance` to reduce the time, as:

\(A, [zeros(∂(d)); f]; tolerance=1E-6);

include("PDE_Helmholtz.jl")

# Plot the solution
import Plots
Plots.plot(u; st=:surface, legend=false)
