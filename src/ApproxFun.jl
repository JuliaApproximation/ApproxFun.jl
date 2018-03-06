__precompile__()

module ApproxFun
    using Base, RecipesBase, FastGaussQuadrature, FastTransforms, DualNumbers,
            BlockArrays, BandedMatrices, BlockBandedMatrices, IntervalSets, Compat
    import StaticArrays, ToeplitzMatrices, Calculus

if VERSION < v"0.7-"
    using Base.FFTW
    import Base.FFTW: Plan
    import Base.FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                    r2r!, r2r
else
    using AbstractFFTs, FFTW
    import AbstractFFTs: Plan
    import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                    r2r!, r2r
end

import Base.LinAlg: BlasInt, BlasFloat, norm, A_ldiv_B!

import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, start, next, done,
                >=, /, ^, \, ∪, transpose, size, to_indexes, reindex, tail, broadcast, broadcast!, copy!, to_index


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


import BlockArrays: nblocks, blocksize, global2blockindex, globalrange, BlockSizes

import BandedMatrices: bandinds, bandrange, PrintShow, bandshift,
                        inbands_getindex, inbands_setindex!, bandwidth, AbstractBandedMatrix,
                        dot, dotu, normalize!, flipsign,
                        colstart, colstop, colrange, rowstart, rowstop, rowrange,
                        bandwidths, showarray, _BandedMatrix, BandedMatrix

import BlockBandedMatrices: blockbandwidth, blockbandwidths, blockcolstop, blockcolrange,
                            blockcolstart, blockrowstop, blockrowstart, blockrowrange,
                            subblockbandwidth, subblockbandwidths, _BlockBandedMatrix,
                            _BandedBlockBandedMatrix, BandedBlockBandedMatrix, BlockBandedMatrix,
                            isblockbanded, isbandedblockbanded, bb_numentries, BlockBandedSizes,
                            BandedBlockBandedSizes

# convenience for 1-d block ranges
const BlockRange1 = BlockRange{1,Tuple{UnitRange{Int}}}

import Base: view

import StaticArrays: SVector

import AbstractFFTs: Plan

if VERSION < v"0.7-"
    using Base.FFTW
    import Base.FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                    r2r!, r2r
else
    using FFTW
    import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                    r2r!, r2r
end


const Vec{d,T} = SVector{d,T}

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
