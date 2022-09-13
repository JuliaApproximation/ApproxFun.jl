# # Solving PDEs

# ## Poisson equation

# We solve the Poisson equation
# ```math
# Δu(x,y) = 0
# ```
# on the rectangle `-1..1 × -1..1`,
# with zero boundary conditions:

using ApproxFun
using LinearAlgebra

d = (-1..1)^2
x,y = Fun(d)
f = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D
A = [Dirichlet(d); Laplacian()]
u = A \ [zeros(∂(d)); f];

# Using a QR Factorization reduces the cost of subsequent calls substantially

QR = qr(A)
u = QR \ [zeros(∂(d)); f];

import Plots
Plots.plot(u; st=:surface, legend=false)

# Many PDEs have weak singularities at the corners,
# in which case it is beneficial to specify a `tolerance` to reduce the time, as:

\(A, [zeros(∂(d)); f]; tolerance=1E-6);

# ## Helmholtz Equation

# We solve
# ```math
# (Δ + 100)u(x,y) = 0
# ```
# on the rectangle `-1..1 × -1..1`,
# subject to the condition that ``u=1`` on the boundary of the domain.

include("PDE1.jl")

# Plot the solution
import Plots
Plots.plot(u; st=:surface, legend=false)
