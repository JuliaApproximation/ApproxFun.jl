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
