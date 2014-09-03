if isdir(Pkg.dir("Gadfly"))
    include("Gadfly.jl")
end
if isdir(Pkg.dir("GLPlot"))
    include("GLPlot.jl")
end
