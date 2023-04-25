# Construct the domain
T = 12
d = 0..T;

Dt = Derivative(d);
ζ = 0.2 # damping ratio
ω0 = 2 # oscillation frequency
L = Dt^2 + 2ζ * ω0 * Dt + ω0^2;

# initial conditions
y0 = 4; # initial displacement
dty0 = 3; # initial velocity

# The differential operator along with the initial condition evaluation operator
A = [L; ivp()];

# We solve the differential equation
y = \(A, [0,y0,dty0], tolerance=1e-6);
