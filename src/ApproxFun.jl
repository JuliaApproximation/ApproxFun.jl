module ApproxFun
    using Base

export Fun,IFun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate

export domain

import Base.values



include("LinearAlgebra/LinearAlgebra.jl")

include("Domains/Domains.jl")
include("Spaces/FunctionSpace.jl")
include("Fun/Fun.jl")
include("Fun/constructors.jl")


include("Operators/Operator.jl")
include("Multivariate/Multivariate.jl")


include("Spaces/Spaces.jl")




## Further extra features

include("Singularities/Singularities.jl")
include("PDE/PDE.jl")
include("Plot/Plot.jl")
include("Extras/Extras.jl")










end #module


