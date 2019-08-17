module ApproxFun
using Base, Reexport, BlockArrays, BandedMatrices, BlockBandedMatrices, DomainSets, IntervalSets,
		SpecialFunctions, AbstractFFTs, FFTW, SpecialFunctions, DSP, DualNumbers, FastTransforms,
		LinearAlgebra, SparseArrays, LowRankApprox, FillArrays, InfiniteArrays, RecipesBase #, Arpack

import Calculus

@reexport using ApproxFunBase    
@reexport using ApproxFunFourier
@reexport using ApproxFunOrthogonalPolynomials
@reexport using ApproxFunSingularities

import ApproxFunBase: normalize!, flipsign, FiniteRange, Fun, MatrixFun, UnsetSpace, VFun, RowVector,
                    UnivariateSpace, AmbiguousSpace, SumSpace, SubSpace, WeightSpace, NoSpace, Space,
                    HeavisideSpace, PointSpace,
                    IntervalOrSegment, RaggedMatrix, AlmostBandedMatrix,
                    AnyDomain, ZeroSpace, ArraySpace, TrivialInterlacer, BlockInterlacer, 
                    AbstractTransformPlan, TransformPlan, ITransformPlan,
                    ConcreteConversion, ConcreteMultiplication, ConcreteDerivative, ConcreteIntegral, CalculusOperator,
                    ConcreteVolterra, Volterra, VolterraWrapper,
                    MultiplicationWrapper, ConversionWrapper, DerivativeWrapper, Evaluation, EvaluationWrapper,
                    Conversion, defaultConversion, defaultcoefficients, default_Fun, Multiplication, Derivative, Integral, bandwidths, 
                    ConcreteEvaluation, ConcreteDefiniteLineIntegral, ConcreteDefiniteIntegral, ConcreteIntegral,
                    DefiniteLineIntegral, DefiniteIntegral, ConcreteDefiniteIntegral, ConcreteDefiniteLineIntegral, IntegralWrapper,
                    ReverseOrientation, ReverseOrientationWrapper, ReverseWrapper, Reverse, NegateEven, 
                    Dirichlet, ConcreteDirichlet, DirichletWrapper,
                    TridiagonalOperator, SubOperator, Space, @containsconstants, spacescompatible,
                    hasfasttransform, canonicalspace, domain, setdomain, prectype, domainscompatible, 
                    plan_transform, plan_itransform, plan_transform!, plan_itransform!, transform, itransform, hasfasttransform, 
                    CanonicalTransformPlan, ICanonicalTransformPlan,
                    Integral, 
                    domainspace, rangespace, boundary, 
                    union_rule, conversion_rule, maxspace_rule, conversion_type, maxspace, hasconversion, points, 
                    rdirichlet, ldirichlet, lneumann, rneumann, ivp, bvp, 
                    linesum, differentiate, integrate, linebilinearform, bilinearform, 
                    UnsetNumber, coefficienttimes, subspace_coefficients, sumspacecoefficients, specialfunctionnormalizationpoint,
                    Segment, IntervalOrSegmentDomain, PiecewiseSegment, isambiguous, Vec, eps, isperiodic,
                    arclength, complexlength,
                    invfromcanonicalD, fromcanonical, tocanonical, fromcanonicalD, tocanonicalD, canonicaldomain, setcanonicaldomain, mappoint,
                    reverseorientation, checkpoints, evaluate, mul_coefficients, coefficients, coefficientmatrix, isconvertible,
                    clenshaw, ClenshawPlan, sineshaw,
                    toeplitz_getindex, toeplitz_axpy!, sym_toeplitz_axpy!, hankel_axpy!, ToeplitzOperator, SymToeplitzOperator, hankel_getindex, 
                    SpaceOperator, ZeroOperator, InterlaceOperator,
                    interlace!, reverseeven!, negateeven!, cfstype, pad!, alternatesign!, mobius,
                    extremal_args, hesseneigvals, chebyshev_clenshaw, recA, recB, recC, roots,splitatroots,
                    chebmult_getindex, intpow, alternatingsum,
                    domaintype, diagindshift, rangetype, weight, isapproxinteger, default_Dirichlet, scal!, dotu,
                    components, promoterangespace, promotedomainspace, choosedomainspace,
                    block, blockstart, blockstop, blocklengths, isblockbanded, pointscompatible,
                    AbstractProductSpace, MultivariateFun, BivariateSpace, 
                    @wrapperstructure, @wrapperspaces, @wrapper, @calculus_operator, resizedata!, slnorm,
                    sample, chop!

import ApproxFunOrthogonalPolynomials: order

import DomainSets: Domain, indomain, UnionDomain, ProductDomain, FullSpace, Point, elements, DifferenceDomain,
            Interval, ChebyshevInterval, boundary, ∂, rightendpoint, leftendpoint,
            dimension, EuclideanDomain

import AbstractFFTs: Plan, fft, ifft
import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                r2r!, r2r,  plan_fft, plan_ifft, plan_ifft!, plan_fft!

import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, iterate,
                >=, /, ^, \, ∪, transpose, size, tail, broadcast, broadcast!, copyto!, copy, to_index, (:),
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



import Base: view

import StaticArrays: StaticArray, SVector


import IntervalSets: (..), endpoints


##Testing
export bisectioninv


## Further extra features

include("Extras/Extras.jl")
include("Plot/Plot.jl")
include("docs.jl")

#include("precompile.jl")
#_precompile_()

end #module
