# # System of equations

# Systems of equations can be handled by creating a matrix of operators and functionals.
# For example, we can solve the system

# ```math
# \begin{gathered}
#     \mathop{u}'' - \mathop{u} + 2 \mathop{v} = \mathop{e}^x,  \\
#     \mathop{u} + \mathop{v}' + \mathop{v} = \cos{x}, \\
#     \mathop{u}(-1) = \mathop{u}'(-1) = \mathop{v}(-1) = 0
# \end{gathered}
# ```

# using the following code:

include("System1.jl")

import Plots
Plots.plot(u, label="u", xlabel="x", legend=:topleft)
Plots.plot!(v, label="v")

# In this example, the automatic space detection failed and so we needed
# to specify explicitly that the domain space for `B` is `Chebyshev()`.
