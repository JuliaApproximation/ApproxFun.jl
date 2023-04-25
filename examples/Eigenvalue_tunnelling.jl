# ### Tunnelling

# We solve the Schrodinger equation in an infinite square well in the domain `-Lx/2..Lx/2`, with a barrier at the center from `-d/2..d/2`
# ```math
# \mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + V(x),
# ```
# where we choose a smoothed barrier potential ``V(x)``.

Lx = 4 # size of the domain
d = 1 # size of the barrier
Δ = 0.1 # smoothing scale of the barrier
x = Fun(-Lx/2..Lx/2)
V = 5 * (1 + tanh((x + d/2)/Δ)) * (1 - tanh((x - d/2)/Δ));
H = -𝒟^2/2 + V;

S = space(x)
B = Dirichlet(S);

λ, v = ApproxFun.eigs(B, H, 1000, tolerance=1E-8);
