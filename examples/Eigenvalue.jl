# # Eigenvalue problem

# ## Standard eigenvalue problem

# In analogy to linear algebra, many differential equations may be posed as eigenvalue problems.
# That is, for some differential operator ``\mathop{L}``, there are a family of functions
# ``\mathop{u}_i(x)`` such that
# ```math
# \mathop{L} \mathop{u}_i(x) = λ_i \mathop{u}_i(x),
# ```
# where ``λ_i`` is the ``i^{th}`` eigenvalue of the ``L`` and has a corresponding
# *eigenfunction* ``\mathop{u}_i(x)``.
# A classic eigenvalue problem is known as the quantum harmonic oscillator where
# ```math
# \mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + \frac{1}{2} x^2,
# ```
# and one demands that ``\mathop{u}(∞) = \mathop{u}(-∞) = 0``.
# Because we expect the solutions to be exponentially suppressed for large ``x``,
# we can approximate this with Dirichlet boundary conditions at a 'reasonably large' ``x``
# without much difference.

# We can express this in ApproxFun as the following:
using ApproxFun
using LinearAlgebra

x = Fun(-8..8)
V = x^2/2
L = -𝒟^2/2 + V
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 500,tolerance=1E-10);

# Note that boundary conditions must be specified in the call to `eigs`.
# Plotting the first ``20`` eigenfunctions offset vertically by their eigenvalue, we see

import Plots
p = Plots.plot(V; legend=false, ylim=(-Inf, λ[22]))
for k=1:20
    Plots.plot!(real(v[k]/norm(v[k]) + λ[k]))
end
p

# If the solutions are not relatively constant near the boundary then one should push
# the boundaries further out.

# For problems with different contraints or boundary conditions,
# `B` can be any zero functional constraint, e.g., `DefiniteIntegral()`.

# ## Self-adjoint Eigenvalue Problem
# Ref:
# [J. L. Aurentz & R. M. Slevinsky (2019), arXiv:1903.08538](https://arxiv.org/abs/1903.08538)

# We solve the confined anharmonic oscillator
# ```math
# \left[-\frac{d^2}{dx^2} + V(x)\right] u = λu,
# ```
# where ``u(\pm 10) = 0``, ``V(x) = ωx^2 + x^4``, and ``ω = 25``.

# Define parameters
ω = 25.0
d = -10..10;
S = Legendre(d) # Equivalently, Ultraspherical(0.5, d)
NS = NormalizedPolynomialSpace(S) # NormalizedLegendre
D1 = Conversion(S, NS)
D2 = Conversion(NS, S);

# Construct the differential operator
V = Fun(x -> ω * x^2 + x^4, S)
L = -Derivative(S, 2) + V;

# Basis recombination to transform to one that satisfies Dirichlet boundary conditions
B = Dirichlet(S)
QS = QuotientSpace(B)
Q = Conversion(QS, S)
R = D1*Q;

# This inversion is computed approximately, such that
# ``\mathrm{C}^{-1} \mathrm{C} ≈ \mathrm{I}`` up to a certain bandwidth
C = Conversion(domainspace(L), rangespace(L))
P = cache(PartialInverseOperator(C, (0, ApproxFun.bandwidth(L, 1) + ApproxFun.bandwidth(R, 1) + ApproxFun.bandwidth(C, 2))));

A = R'D1*P*L*D2*R
B = R'R;

# We impose a cutoff to obtain approximate finite matrix representations
n = 3000
SA = Symmetric(A[1:n,1:n], :L)
SB = Symmetric(B[1:n,1:n], :L)
λ = eigvals(SA, SB)[1:round(Int, 3n/5)];

# We plot the eigenvalues
import Plots
Plots.plot(λ; title = "Eigenvalues", legend=false)
