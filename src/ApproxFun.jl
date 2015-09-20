VERSION >= v"0.4.0-dev+6641" && __precompile__()

module ApproxFun
    using Base, Compat


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
