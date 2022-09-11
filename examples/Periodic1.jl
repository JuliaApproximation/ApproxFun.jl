using ApproxFun
using LinearAlgebra

# We may approximate a periodic function in `Fourier` space
f = Fun(cos, Fourier(-π..π))
norm(differentiate(f) + Fun(sin, Fourier(-π..π))) < 100eps()

# ## Boundary Value Problems

# Due to the periodicity, Fourier representations allow for the asymptotic
# savings of 2/π in the number of coefficients that need
# to be stored compared with a Chebyshev representation.

# ODEs can also be solved when the solution is periodic.

# In this example, we solve
# ```math
# \frac{dy}{dt} + (1 + \sin(\cos(2t)))y = \exp(\sin(10t))
# ```
# subject to periodic boundary conditions on the domain ``[-\pi,\,\pi]``.

# Solve the differential equation in Chebyshev space
s = Chebyshev(-π..π)
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
B = periodic(s,0)
uChebyshev = [B;L]\[0.,f]
ncoefficients(uChebyshev)

# Solve the differential equation in Fourier space
s = Fourier(-π..π)
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
uFourier = L\f
ncoefficients(uFourier)

using Test #src
@test uChebyshev(0.) ≈ uFourier(0.) #src
