# # Sampling
include("Sampling1.jl")

# We plot a histogram of the samples
import Plots
Plots.histogram(x; bins=100, normalize = :pdf,
    title = "Histogram", xlabel="value", ylabel="frequency",
    xlim = (-5,5), ylim=(0, Inf), legend = false)
Plots.plot!(f/sum(f), linewidth=3)
