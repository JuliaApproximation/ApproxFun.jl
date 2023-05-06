# # Eigenvalue problem

include("Eigenvalue_standard.jl")

# Note that boundary conditions must be specified in the call to `eigs`.
# Plotting the first ``20`` eigenfunctions offset vertically by their eigenvalue, we see

import Plots
using LinearAlgebra: norm
p = Plots.plot(V, legend=false, ylim=(-Inf, λ[22]))
for k=1:20
    Plots.plot!(real(v[k]) + real(λ[k]))
end
p

# If the solutions are not relatively constant near the boundary then one should push
# the boundaries further out.

# For problems with different contraints or boundary conditions,
# `B` can be any zero functional constraint, e.g., `DefiniteIntegral()`.

include("Eigenvalue_anharmonic.jl")

# We plot the eigenvalues
import Plots
Plots.plot(λ, title = "Eigenvalues", legend=false)

include("Eigenvalue_well_barrier.jl")

# We plot the first few eigenfunctions offset by their eigenvalues.
# The eigenfunctions appear in odd-even pairs by construction.

import Plots
using LinearAlgebra: norm
p1 = Plots.plot(Vfull, legend=false, ylim=(-Inf, λ[14]), title="Eigenfunctions",
    seriestype=:path, linestyle=:dash, linewidth=2)
Plots.vline!([-Lx/2, Lx/2], color=:black)
for k=1:12
    Plots.plot!(real(v[k]) + λ[k])
end

# Zoom into the ground state:
p2 = Plots.plot(Vfull, legend=false, ylim=(-Inf, λ[3]), title="Ground state",
    seriestype=:path, linestyle=:dash, linewidth=2)
Plots.vline!([-Lx/2, Lx/2], color=:black)
Plots.plot!(real(v[1]) + λ[1], linewidth=2)

Plots.plot(p1, p2, layout = 2)
