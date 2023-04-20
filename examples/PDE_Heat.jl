# ## Heat equation

# We solve the heat equation
# ```math
# \frac{\partial}{\partial t} u(x,t) = \frac{\partial^2}{\partial x^2} u(x,t)
# ```
# on the spatial domain `-1..1` from ``t=0`` to ``t=1``,
# with zero boundary conditions, and the initial condition
# ``u(x,0) = u_0(x)``, where we choose ``u_0(x)`` to be a sharp Gaussian.

using ApproxFun
using LinearAlgebra

dx = -1..1;
dt = 0..1;

# We construct a 2D domain that is the tensor product of x and t
d = dx × dt;

# The initial condition
u0 = Fun(x->exp(-x^2/(2*0.2^2)), dx);

# We define the derivatives on the 2D domain
Dx = Derivative(d,[1,0]);
Dt = Derivative(d,[0,1]);

# The operator along with the initial and the boundary conditions
A = [Dt - Dx^2; I⊗ldirichlet(dt); bvp(dx)⊗I];

# We collect the initial condition along with zeros for the two boundary conditions and the equation
rhs = [0.0; u0; 0.0; 0.0];

# Solve the equation with an arbitrary tolerance.
# A lower tolerance will increase accuracy at the expense of execution time
u = \(A, rhs, tolerance=1e-4);
