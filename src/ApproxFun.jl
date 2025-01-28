__precompile__()

module ApproxFun
    using Base, Compat, Plots, FastGaussQuadrature, FastTransforms, DualNumbers, BandedMatrices
    import FixedSizeArrays, ToeplitzMatrices

import Base.LinAlg: BlasFloat

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
                >=,./,/,.^,^,\,∪,transpose, size


import BandedMatrices: bzeros, bandinds, bandrange, PrintShow, eachbandedindex, bandshift,
                        unsafe_getindex, unsafe_setindex!, bandwidth, AbstractBandedMatrix,
                        dot, dotu, normalize!, flipsign,
                        colstart, colstop, colrange, rowstart, rowstop, rowrange,
                        bandwidths, αA_mul_B_plus_βC!, showarray

import Compat: view

import FixedSizeArrays: Vec

export pad!,pad,sample,chop!,complexroots,roots,svfft, reverseorientation

##Testing
export bisectioninv





include("LinearAlgebra/LinearAlgebra.jl")


include("Fun/Fun.jl")
include("Multivariate/Multivariate.jl")
include("Operators/Operator.jl")

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")




## Further extra features

include("PDE/PDE.jl")
include("Caching/caching.jl")
include("Extras/Extras.jl")
include("Plot/Plot.jl")

include("docs.jl")



## precompile
function _precompile_()
    precompile(Chebyshev,tuple())
    precompile(Fun,tuple())
    precompile(Base.exp,(Fun{Chebyshev{Interval{Float64}},Float64},))
    precompile(linsolve,(Vector{Operator{Float64}},Vector{Float64}))
    precompile(\,(Vector{Operator{Float64}},Vector{Float64}))
    precompile(+,(Int,Fun{Chebyshev{Interval{Float64}},Float64}))
    precompile(+,(Fun{Chebyshev{Interval{Float64}},Float64},Fun{Chebyshev{Interval{Float64}},Float64}))
    precompile(+,(Fun{ConstantSpace{AnyDomain},Float64},Fun{Chebyshev{Interval{Float64}},Float64}))
    precompile(union,(ConstantSpace{AnyDomain},Chebyshev{Interval{Float64}}))
    precompile(union,(Chebyshev{Interval{Float64}},ConstantSpace{AnyDomain}))
    precompile(Fun,(Fun{ConstantSpace{AnyDomain},Float64},Chebyshev{Interval{Float64}}))
end

_precompile_()

end #module
