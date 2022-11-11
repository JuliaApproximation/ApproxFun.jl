module ApproxFun
using Base, Reexport,
		AbstractFFTs, FFTW, DualNumbers, FastTransforms,
		LinearAlgebra, RecipesBase, DomainSets, SpecialFunctions #, Arpack

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
                    sample, chop!, isbanded, colrange, bandwidth

import ApproxFunOrthogonalPolynomials: order


import AbstractFFTs: Plan, fft, ifft
import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                r2r!, r2r,  plan_fft, plan_ifft, plan_ifft!, plan_fft!

import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, iterate,
                >=, /, ^, \, ∪, transpose, size, tail, broadcast, broadcast!, copyto!, copy, to_index, (:),
                similar, map, vcat, hcat, hvcat, show, summary, stride, sum, cumsum, sign, conj, inv,
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

# import Arpack: eigs


import FastTransforms: ChebyshevTransformPlan, IChebyshevTransformPlan, plan_chebyshevtransform,
                        plan_chebyshevtransform!, plan_ichebyshevtransform, plan_ichebyshevtransform!

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


##Testing
export bisectioninv


## Further extra features

include("Extras/Extras.jl")
include("Plot/Plot.jl")

#include("precompile.jl")
#_precompile_()

end #module
