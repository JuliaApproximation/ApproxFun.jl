__precompile__()

module ApproxFun
    using Base, Compat, Plots, FastGaussQuadrature, FastTransforms, DualNumbers, BandedMatrices
    import FixedSizeArrays, ToeplitzMatrices, Calculus

import Base.LinAlg: BlasFloat

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,!,!=,
                >=,./,/,.^,^,\,∪,transpose, size

# we need to import all special functions to use Calculus.symbolic_derivatives_1arg
# we can't do importall Base as we replace some Base definitions
import Base: sinpi, cospi, airy, besselh, exp,
                    asinh, acosh,atanh, erfcx, dawson, erf, erfi,
                    sin, cos, sinh, cosh, airyai, airybi, airyaiprime, airybiprime,
                    hankelh1, hankelh2, besselj, bessely, besseli, besselk,
                    besselkx, hankelh1x, hankelh2x, exp2, exp10, log2, log10,
                    tan, tanh, csc, asin, acsc, sec, acos, asec,
                    cot, atan, acot, sinh, csch, asinh, acsch,
                    sech, acosh, asech, tanh, coth, atanh, acoth,
                    expm1, log1p, lfact, sinc, cosc, erfinv, erfcinv, beta, lbeta,
                    eta, zeta, gamma,  lgamma, polygamma, invdigamma, digamma, trigamma,
                    abs, sign, log, expm1, tan, abs2, sqrt, angle, max, min, cbrt, log,
                    atan, acos, asin, erfc, inv


import BandedMatrices: bzeros, bandinds, bandrange, PrintShow, eachbandedindex, bandshift,
                        unsafe_getindex, unsafe_setindex!, bandwidth, AbstractBandedMatrix,
                        dot, dotu, normalize!, flipsign,
                        colstart, colstop, colrange, rowstart, rowstop, rowrange,
                        bandwidths, αA_mul_B_plus_βC!, showarray

import Compat: view

import FixedSizeArrays: Vec

export pad!, pad, chop!, sample,
       complexroots, roots, svfft, isvfft,
       reverseorientation

##Testing
export bisectioninv





include("LinearAlgebra/LinearAlgebra.jl")


include("Fun/Fun.jl")

include("Domains/Domains.jl")

include("Multivariate/Multivariate.jl")
include("Operators/Operator.jl")

include("Spaces/Spaces.jl")




## Further extra features

include("PDE/PDE.jl")
include("Caching/caching.jl")
include("Extras/Extras.jl")
include("Plot/Plot.jl")

include("docs.jl")

include("testing.jl")


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
