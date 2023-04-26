# ## Self-adjoint Eigenvalue Problem
# Ref:
# [J. L. Aurentz & R. M. Slevinsky (2019), arXiv:1903.08538](https://arxiv.org/abs/1903.08538)

# We solve the confined anharmonic oscillator
# ```math
# \left[-\frac{d^2}{dx^2} + V(x)\right] u = λu,
# ```
# where ``u(\pm 10) = 0``, ``V(x) = ωx^2 + x^4``, and ``ω = 25``.

using ApproxFun
using LinearAlgebra
using BandedMatrices

# Define parameters
ω = 25.0
d = -10..10;
S = Legendre(d); # Equivalently, Ultraspherical(0.5, d)

# Construct the differential operator
V = Fun(x -> ω * x^2 + x^4, S)
L = -Derivative(S, 2) + V;

# Boundary conditions that are used in the basis recombination
B = Dirichlet(S);

# The system may be recast as a generalized eigenvalue problem
# ``A_S\,v = λ\, B_S\, v``, where ``A_S`` and ``B_S`` are symmetric band matrices.
# We wrap the operators in `SymmetricEigensystem` to implicitly perform the basis recombination

SEg = ApproxFun.SymmetricEigensystem(L, B);

# We construct `n × n` matrix representations of the opertors that we diagonalize
n = 3000
λ = eigvals(SEg, n);

# We retain a fraction of the eigenvalues with the least magnitude
λ = λ[1:round(Int, 3n/5)];
