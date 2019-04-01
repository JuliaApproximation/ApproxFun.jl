__precompile__()

module ApproxFun
    using Base, Reexport, BlockArrays, BandedMatrices, BlockBandedMatrices, DomainSets, IntervalSets,
            SpecialFunctions, AbstractFFTs, FFTW, SpecialFunctions, DSP, DualNumbers, FastTransforms,
            LinearAlgebra, SparseArrays, LowRankApprox, FillArrays, InfiniteArrays #, Arpack

@reexport using ApproxFunBase    
@reexport using ApproxFunFourier
@reexport using ApproxFunOrthogonalPolynomials

import ApproxFunBase: normalize!, flipsign, FiniteRange, MatrixFun, UnsetSpace, VFun, RowVector,
                    UnivariateSpace, AmbiguousSpace, IntervalOrSegment, RaggedMatrix, AlmostBandedMatrix,
                    AnyDomain, ZeroSpace, TrivialInterlacer, BlockInterlacer, TransformPlan, ITransformPlan,
                    ConcreteConversion, ConcreteMultiplication, ConcreteDerivative, ConcreteEvaluation,
                    TridiagonalOperator, SubOperator, Space

import DomainSets: Domain, indomain, UnionDomain, ProductDomain, FullSpace, Point, elements, DifferenceDomain,
            Interval, ChebyshevInterval, boundary, ∂, rightendpoint, leftendpoint,
            dimension, Domain1d, Domain2d




import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, iterate,
                >=, /, ^, \, ∪, transpose, size, reindex, tail, broadcast, broadcast!, copyto!, copy, to_index, (:),
                similar, map, vcat, hcat, hvcat, show, summary, stride, sum, cumsum, sign, imag, conj, inv,
                complex, reverse, exp, sqrt, abs, abs2, sign, issubset, values, in, first, last, rand, intersect, setdiff,
                isless, union, angle, join, isnan, isapprox, isempty, sort, merge, promote_rule,
                minimum, maximum, extrema, argmax, argmin, findmax, findmin, isfinite,
                zeros, zero, one, promote_rule, repeat, length, resize!, isinf,
                getproperty, findfirst, unsafe_getindex, fld, cld, div, real, imag,
                @_inline_meta, eachindex, firstindex, lastindex, keys, isreal, OneTo,
                Array, Vector, Matrix, view, ones, @propagate_inbounds, print_array,
                split

import Base.Broadcast: BroadcastStyle, Broadcasted, AbstractArrayStyle, broadcastable,
                        DefaultArrayStyle, broadcasted



import LinearAlgebra: BlasInt, BlasFloat, norm, ldiv!, mul!, det, eigvals, dot, cross,
                        qr, qr!, rank, isdiag, istril, istriu, issymmetric, ishermitian,
                        Tridiagonal, diagm, diagm_container, factorize, nullspace,
                        Hermitian, Symmetric, adjoint, transpose, char_uplo

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

import StaticArrays: SVector

import BlockArrays: nblocks, blocksize, global2blockindex, globalrange, BlockSizes

import BandedMatrices: bandrange, bandshift,
                        inbands_getindex, inbands_setindex!, bandwidth, AbstractBandedMatrix,
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
import LazyArrays: cache
import InfiniteArrays: Infinity, InfRanges, AbstractInfUnitRange, OneToInf

"""
`Curve` Represents a domain defined by the image of a Fun.  Example
usage would be

```julia
x=Fun(1..2)
Curve(exp(im*x))  # represents an arc
```
"""
const Curve{S,T} = Union{IntervalCurve{S,T},PeriodicCurve{S,T}}
Curve(f::Fun{<:Space{<:PeriodicDomain}}) = PeriodicCurve(f)

#TODO: Make type stable
Curve(f::Fun{<:Space{<:ChebyshevInterval}}) = IntervalCurve(f) 

export Curve

const AffineDomain = Union{AbstractInterval,Segment,PeriodicSegment,Ray,Line}

## set minus
function Base.setdiff(d::AffineDomain,ptsin::UnionDomain{AS}) where {AS <: AbstractVector{P}} where {P <: Point}
    pts=Number.(elements(ptsin))
    isempty(pts) && return d
    tol=sqrt(eps(arclength(d)))
    da=leftendpoint(d)
    isapprox(da,pts[1];atol=tol) && popfirst!(pts)
    isempty(pts) && return d
    db=rightendpoint(d)
    isapprox(db,pts[end];atol=tol) && pop!(pts)

    sort!(pts)
    leftendpoint(d) > rightendpoint(d) && reverse!(pts)
    filter!(p->p ∈ d,pts)

    isempty(pts) && return d
    length(pts) == 1 && return d \ pts[1]

    ret = Array{Domain}(undef, length(pts)+1)
    ret[1] = Domain(leftendpoint(d) .. pts[1])
    for k = 2:length(pts)
        ret[k] = Domain(pts[k-1]..pts[k])
    end
    ret[end] = Domain(pts[end] .. rightendpoint(d))
    UnionDomain(ret)
end

# convenience for 1-d block ranges
const BlockRange1 = BlockRange{1,Tuple{UnitRange{Int}}}

import Base: view

import StaticArrays: StaticArray, SVector


import IntervalSets: (..), endpoints

const Vec{d,T} = SVector{d,T}

export pad!, pad, chop!, sample,
       complexroots, roots, svfft, isvfft,
       reverseorientation, jumplocations

##Testing
export bisectioninv

export .., Interval, ChebyshevInterval, leftendpoint, rightendpoint, endpoints



include("Spaces/Spaces.jl")




## Further extra features

include("Extras/Extras.jl")
include("Plot/Plot.jl")
include("docs.jl")

#include("precompile.jl")
#_precompile_()

end #module
