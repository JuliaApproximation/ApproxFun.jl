using ApproxFun
using LinearAlgebra

# The rectangular domain may be expressed as a tensor product of `ChebyshevInterval`s
d = ChebyshevInterval()^2;

# We construct the operator that collates the differential operator and boundary conditions
L = [Dirichlet(d); Laplacian()+100I];

# We compute the QR decomposition of the operator to speed up the solution
Q = qr(L);
ApproxFun.resizedata!(Q,:,4000);

# The boundary condition is a function that is equal to one on each edge
boundary_cond = ones(∂(d));

# Solve the differential equation.
# This function has weak singularities at the corner,
# so we specify a lower tolerance to avoid resolving these singularities completely
u = \(Q, [boundary_cond; 0.]; tolerance=1E-7);

#src # We verify that the solution satisfies the boundary condition
using Test #src
@test u(0.1,1.) ≈ 1.0 #src

#src # We validate the solution at an internal point
@test ≈(u(0.1,0.2), -0.02768276827514463; atol=1E-8) #src
