module ApproxFun
    using Base, Compat

export Fun,IFun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate

export domain,space

import Base.values



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










end #module


