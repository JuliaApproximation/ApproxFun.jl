# # Self-adjoint Eigenvalue Problem
# Ref:
# [J. L. Aurentz & R. M. Slevinsky (2019), arXiv:1903.08538](https://arxiv.org/abs/1903.08538)

# We solve the confined anharmonic oscillator
# ```math
# \left[-\frac{d^2}{dx^2} + V(x)\right] u = λu,
# ```
# where ``u(\pm 10) = 0``, ``V(x) = ωx^2 + x^4``, and ``ω = 25``.

using ApproxFun
using LinearAlgebra

# Define parameters
ω = 25.0
d = -10..10;
S = Legendre(d)
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
# $\mathrm{C}^-1 * \mathrm{C} ≈ \mathrm{I}$ up to a certain bandwidth
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
