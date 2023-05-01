# ### Infinite well with a barrier

# We solve the Schrodinger equation in an infinite square well in the domain `-Lx/2..Lx/2`, with a finite barrier in the middle from `-d/2..d/2`
# ```math
# \mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + V(x),
# ```
# where ``V(x)`` is given by an infinite well with a smoothed barrier at the center.

# We note that the system has parity symmetry, so the solutions may be separated into odd and even functions.
# We may therefore solve the problem only for half the domain, with Dirichlet boundary conditions at the midpoint
# for odd functions, and Neumann conditions for even functions.
# This projection to subspaces allows us to halve the sizes of matrices that we need to diagonalize

Lx = 4 # size of the domain
S = Legendre(0..Lx/2)
x = Fun(S)
d = 1 # size of the barrier
Δ = 0.1 # smoothing scale of the barrier
V = 50 * (1 - tanh((x - d/2)/Δ))/2; # right half of the potential barrier
H = -Derivative(S)^2/2 + V;

# Odd solutions, with a zero Dirichlet condition at `0` representing a node
B = Dirichlet(S);
Seig = ApproxFun.SymmetricEigensystem(H, B);

# Diagonalize `n × n` matrix representations of the basis-recombined operators
# We specify a tolerance to reject spurious solutions arising from the discretization
n = 100
λodd, vodd = ApproxFun.eigs(Seig, n, tolerance=1e-8);

# To extend the solutions to the full domain, we construct the left-half space.
Scomplement = Legendre(-Lx/2..0);

# We use the fact that Legendre polynomials of odd orders are odd functions,
# and those of even orders are even functions
# Using this, for the odd solutions, we negate the even-order coefficients to construct the odd image in `-Lx/2..0`
oddimage(f, Scomplement) = Fun(Scomplement, [(-1)^isodd(m) * c for (m,c) in enumerate(coefficients(f))]);
voddimage = oddimage.(vodd, Scomplement);

# construct the functions over the entire domain `-Lx/2..Lx/2` as piecewise sums over the two half domains `-Lx/2..0` and `0..Lx/2`
voddfull = voddimage .+ vodd;

# Even solutions, with a Neumann condition at `0` representing the symmetry of the function
B = [lneumann(S); rdirichlet(S)];

Seig = ApproxFun.SymmetricEigensystem(H, B);
λeven, veven = ApproxFun.eigs(Seig, n, tolerance=1e-8);

# For the even solutions, we negate the odd-order coefficients to construct the even image in `-Lx/2..0`
evenimage(f, Scomplement) = Fun(Scomplement, [(-1)^iseven(m) * c for (m,c) in enumerate(coefficients(f))]);
vevenimage = evenimage.(veven, Scomplement);
vevenfull = vevenimage .+ veven;

# We interlace the eigenvalues and eigenvectors to obtain the entire spectrum
λ = vec(permutedims([λeven λodd]));
v = vec(permutedims([vevenfull voddfull]));

# Construct the potential on the full domain for plotting using an even image
Vevenimage = evenimage(V, Scomplement);
Vfull = Vevenimage + V;
