# # Sampling
# Random number sampling using [Olver & Townsend 2013](https://arxiv.org/abs/1307.1223).
# The following code samples 10,000 from a PDF given as the absolute value of the sine function on [-5,5]:
using ApproxFun

f = abs(Fun(sin,-5..5))
x = ApproxFun.sample(f, 10000);

# We plot a histogram of the samples
import Plots
Plots.histogram(x; bins=100, normalize = :pdf,
    title = "Histogram", xlabel="value", ylabel="frequency",
    xlim = (-5,5), ylim=(0, Inf), legend = false)
