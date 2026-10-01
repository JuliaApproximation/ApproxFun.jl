## This includes extra tests that are too time consuming for Travis


include("runtests.jl")


println("Full PDE tests")

include("FullPDETest.jl")


println("Speed tests")
include("SpeedTest.jl")
include("SpeedODETest.jl")
include("SpeedPDETest.jl")



println("Example tests")
if isdir(Pkg.dir("GR")) || isdir(Pkg.dir("Plotly")) || isdir(Pkg.dir("PlotlyJS")) || isdir(Pkg.dir("PyPlot"))
    include("ExamplesTest.jl")
else
    warn("Unable to do Examples since Gadfly.jl is not installed")
end


println("Readme tests")
include("ReadmeTest.jl")
