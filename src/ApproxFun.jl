__precompile__()

module ApproxFun
    using Base, RecipesBase, FastGaussQuadrature, FastTransforms, DualNumbers, BandedMatrices, IntervalSets
    import FixedSizeArrays, ToeplitzMatrices, Calculus

import Base.LinAlg: BlasInt, BlasFloat

import Base: values, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, start, next, done,
                >=, /, ^, \, ∪, transpose, size, to_indexes, reindex, tail, broadcast

if VERSION < v"0.6.0-dev"
    import Base: .*, .+, .-, .^, ./
end

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


import BandedMatrices: bzeros, bandinds, bandrange, PrintShow, bandshift,
                        inbands_getindex, inbands_setindex!, bandwidth, AbstractBandedMatrix,
                        dot, dotu, normalize!, flipsign,
                        colstart, colstop, colrange, rowstart, rowstop, rowrange,
                        bandwidths, αA_mul_B_plus_βC!, showarray

import Base: view

import FixedSizeArrays: Vec

export pad!, pad, chop!, sample,
       complexroots, roots, svfft, isvfft,
       reverseorientation

##Testing
export bisectioninv

export ..



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

#include("precompile.jl")
#_precompile_()

end #module
