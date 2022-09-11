using ApproxFun
using LinearAlgebra

# Define the vector that collates the differential equation and
# the boundary conditions
N(u, x = Fun()) = [u(-1.)-1., u(1.)+0.5, 0.001u'' + 6(1-x^2)u' + u^2 - 1];

# Solve the equation using Newton iteration
function nbvpsolver()
    x = Fun()
    u0 = 0 * x # starting value

    newton(N, u0)
end

u = nbvpsolver();

#src # We verify that the solution satisfies the differential equation and the boundary conditions
using Test #src
@test norm(N(u)) ≤ 1000eps() #src
