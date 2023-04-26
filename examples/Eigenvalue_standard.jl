# ## Standard eigenvalue problem

# In analogy to linear algebra, many differential equations may be posed as eigenvalue problems.
# That is, for some differential operator ``\mathop{L}``, there are a family of functions
# ``\mathop{u}_i(x)`` such that
# ```math
# \mathop{L} \mathop{u}_i(x) = λ_i \mathop{u}_i(x),
# ```
# where ``λ_i`` is the ``i^{th}`` eigenvalue of the ``L`` and has a corresponding
# *eigenfunction* ``\mathop{u}_i(x)``.

# ### Quantum harmonic oscillator

# A classic eigenvalue problem is known as the quantum harmonic oscillator where
# ```math
# \mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + \frac{1}{2} x^2,
# ```
# and one demands that ``\mathop{u}(∞) = \mathop{u}(-∞) = 0``.
# Because we expect the solutions to be exponentially suppressed for large ``x``,
# we can approximate this with Dirichlet boundary conditions at a 'reasonably large' ``x``
# without much difference.

# We can express this in ApproxFun as the following:
using ApproxFun

x = Fun(-8..8)
V = x^2/2
L = -𝒟^2/2 + V
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 500,tolerance=1E-10);
