__precompile__()

module ApproxFun
    using Base, RecipesBase, FastGaussQuadrature, FastTransforms, DualNumbers,
            BlockArrays, BandedMatrices, BlockBandedMatrices, IntervalSets,
            SpecialFunctions, AbstractFFTs, FFTW, SpecialFunctions,
            LinearAlgebra, LowRankApprox, SparseArrays, FillArrays, InfiniteArrays #, Arpack
    import StaticArrays, ToeplitzMatrices, Calculus


import AbstractFFTs: Plan, fft, ifft
import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                r2r!, r2r,  plan_fft, plan_ifft, plan_ifft!, plan_fft!


import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, iterate,
                >=, /, ^, \, ∪, transpose, size, reindex, tail, broadcast, broadcast!, copyto!, copy, to_index, (:),
                similar, map, vcat, hcat, hvcat, show, summary, stride, sum, cumsum, sign, imag, conj, inv,
                complex, reverse, exp, sqrt, abs, abs2, sign, issubset, values, in, first, last, rand, intersect, setdiff,
                isless, union, angle, join, isnan, isapprox, isempty, sort, merge, promote_rule,
                minimum, maximum, extrema, argmax, argmin, findmax, findmin, isfinite,
                zeros, zero, one, promote_rule, repeat, length, resize!, isinf,
                getproperty, findfirst, unsafe_getindex, fld, cld, div, real, imag,
                @_inline_meta, eachindex, lastindex, keys, isreal,
                Array, Vector, Matrix, view, ones, @propagate_inbounds, print_array

import Base.Broadcast: BroadcastStyle, Broadcasted, AbstractArrayStyle, broadcastable,
                        DefaultArrayStyle, broadcasted


import LinearAlgebra: BlasInt, BlasFloat, norm, ldiv!, mul!, det, eigvals, dot, cross,
                        qr, qr!, isdiag, rank, issymmetric, ishermitian, Tridiagonal,
                        diagm, factorize, nullspace, adjoint, transpose, diagm_container

import SparseArrays: blockdiag

# import Arpack: eigs

# we need to import all special functions to use Calculus.symbolic_derivatives_1arg
# we can't do importall Base as we replace some Base definitions
import SpecialFunctions: sinpi, cospi, airy, besselh,
                    asinh, acosh,atanh, erfcx, dawson, erf, erfi,
                    sin, cos, sinh, cosh, airyai, airybi, airyaiprime, airybiprime,
                    hankelh1, hankelh2, besselj, besselj0, bessely, besseli, besselk,
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
                        dotu, normalize!, flipsign,
                        colstart, colstop, colrange, rowstart, rowstop, rowrange,
                        bandwidths, _BandedMatrix, BandedMatrix

import BlockBandedMatrices: blockbandwidth, blockbandwidths, blockcolstop, blockcolrange,
                            blockcolstart, blockrowstop, blockrowstart, blockrowrange,
                            subblockbandwidth, subblockbandwidths, _BlockBandedMatrix,
                            _BandedBlockBandedMatrix, BandedBlockBandedMatrix, BlockBandedMatrix,
                            isblockbanded, isbandedblockbanded, bb_numentries, BlockBandedSizes,
                            BandedBlockBandedSizes

import FastTransforms: ChebyshevTransformPlan, IChebyshevTransformPlan, plan_chebyshevtransform,
                        plan_chebyshevtransform!, plan_ichebyshevtransform, plan_ichebyshevtransform!

import FillArrays: AbstractFill, getindex_value

import InfiniteArrays: Infinity, InfRanges, AbstractInfUnitRange

# convenience for 1-d block ranges
const BlockRange1 = BlockRange{1,Tuple{UnitRange{Int}}}

import StaticArrays: SVector


const Vec{d,T} = SVector{d,T}

export pad!, pad, chop!, sample,
       complexroots, roots, svfft, isvfft,
       reverseorientation, jumplocations

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
