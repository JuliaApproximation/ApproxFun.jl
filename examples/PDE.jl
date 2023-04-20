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

include("PDE_Heat.jl")

# Plot the solution
import Plots
xplot = -1:0.02:1
p = Plots.plot(xplot, u.(xplot, 0), label="t=0", legend=true, linewidth=2)
for t in [0.05, 0.1, 0.2, 0.5, 0.8]
	Plots.plot!(xplot, u.(xplot, t), label="t=$t")
end
p
