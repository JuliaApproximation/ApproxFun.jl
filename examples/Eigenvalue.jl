# # Eigenvalue problem

include("Eigenvalue_standard.jl")

# Note that boundary conditions must be specified in the call to `eigs`.
# Plotting the first ``20`` eigenfunctions offset vertically by their eigenvalue, we see

import Plots
using LinearAlgebra: norm
p = Plots.plot(V, legend=false, ylim=(-Inf, λ[22]))
for k=1:20
    Plots.plot!(real(v[k]/norm(v[k]) + λ[k]))
end
p

# If the solutions are not relatively constant near the boundary then one should push
# the boundaries further out.

include("Eigenvalue_tunnelling.jl")

# We plot the first few eigenfunctions
p = Plots.plot(V, legend=false)
Plots.vline!([-Lx/2, Lx/2], linecolor=:black)
p_twin = Plots.twinx(p)
for k=1:6
    Plots.plot!(p_twin, real(v[k]/norm(v[k]) + λ[k]), label="$k")
end
p

# Note that the parity symmetry isn't preserved exactly at finite matrix sizes.
# In general, it's better to preserve the symmetry of the operator matrices (see section below),
# and projecting them on to the appropriate subspaces.

# For problems with different contraints or boundary conditions,
# `B` can be any zero functional constraint, e.g., `DefiniteIntegral()`.

include("Eigenvalue_symmetric.jl")

# We plot the eigenvalues
import Plots
Plots.plot(λ; title = "Eigenvalues", legend=false)
