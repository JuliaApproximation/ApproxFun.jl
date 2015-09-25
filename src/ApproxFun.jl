__precompile__()

module ApproxFun
    using Base, FastGaussQuadrature


import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,
                >=,./,/,.^,^,\,∪,transpose


export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,fasttimes

##Testing
export bisectioninv, clenshaw





include("LinearAlgebra/LinearAlgebra.jl")


include("Fun/Fun.jl")
include("Multivariate/Multivariate.jl")
include("Operators/Operator.jl")

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")




## Further extra features

include("PDE/PDE.jl")
include("Plot/Plot.jl")
include("Extras/Extras.jl")

include("docs.jl")








end #module
