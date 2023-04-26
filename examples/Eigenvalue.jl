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

# For problems with different contraints or boundary conditions,
# `B` can be any zero functional constraint, e.g., `DefiniteIntegral()`.

include("Eigenvalue_anharmonic.jl")

# We plot the eigenvalues
import Plots
Plots.plot(λ, title = "Eigenvalues", legend=false)

include("Eigenvalue_tunnelling.jl")

# We plot the first few eigenvalues
Plots.plot(λ, title="Eigenvalues", legend=false)
