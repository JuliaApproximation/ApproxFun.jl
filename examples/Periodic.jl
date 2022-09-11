# # Periodic Domains
include("Periodic1.jl")

# We note that the number of coefficients in the Fourier space is lower than that
# in the Chebyshev space by a factor of approximately π/2
ncoefficients(uFourier)/ncoefficients(uChebyshev),2/π

# We plot the solution
import Plots
Plots.plot(uFourier, xlims=(-pi, pi), legend=false)

