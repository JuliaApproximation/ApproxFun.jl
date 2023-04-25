using ApproxFun
using LinearAlgebra

# Construct the domain
a, b = 0, 2pi
d = a..b;

# We construct the differential operator ``L = d^2/dx^2 + 4``:
D = Derivative(d)
L = D^2 + 4;

# We have not chosen the space explicitly, and the solve chooses it to be `Chebyshev(d)` internally

# We incorporate the initial conditions using `ivp`, which is a shortcut for evaluating the function
# and it's derivative at the left boundary of the domain
A = [L; ivp()];

# Initial conditions
u0 = 0;
dtu0 = 2;

# We solve the differential equation
u = \(A, [0,u0,dtu0], tolerance=1e-6);
