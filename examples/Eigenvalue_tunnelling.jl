# ### Tunnelling

# We solve the Schrodinger equation in an infinite square well in the domain `-Lx/2..Lx/2`, with a barrier at the center from `-d/2..d/2`
# ```math
# \mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + V(x),
# ```
# where we choose a smoothed barrier potential ``V(x)``.

# We note that the system has parity symmetry, so the solutions may be separated into odd and even functions.
# We may therefore solve the problem only for half the domain, with Dirichlet boundary conditions at the midpoint
# for odd functions, and Neumann conditions for even functions.
# This projection to subspaces allows us to halve the sizes of matrices that we need to diagonalize

Lx = 4 # size of the domain
d = 1 # size of the barrier
Δ = 0.1 # smoothing scale of the barrier
S = Legendre(0..Lx/2)
x = Fun(S)
V = 5 * (1 - tanh((x - d/2)/Δ));
H = -Derivative(S)^2/2 + V;

# Odd solutions
B = Dirichlet(S);
SEg = ApproxFun.SymmetricEigensystem(H, B);

# Diagonalize `n × n` matrix representations of the basis-recombined operators
n = 100
λodd = eigvals(SEg, n);
λodd = λodd[1:round(Int, 3n/5)];

# Even solutions
B = [lneumann(S); rdirichlet(S)];

SEg = ApproxFun.SymmetricEigensystem(H, B);
λeven = eigvals(SEg, n);
λeven = λeven[1:round(Int, 3n/5)];

# We interlace the eigenvalues to obtain the entire spectrum
λ = vec([λeven λodd]');
