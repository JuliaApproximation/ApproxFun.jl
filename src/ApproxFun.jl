module ApproxFun
using Base, Reexport,
		AbstractFFTs, FFTW, DualNumbers, FastTransforms,
		LinearAlgebra, RecipesBase, DomainSets, SpecialFunctions #, Arpack

import Calculus

@reexport using ApproxFunBase
@reexport using ApproxFunFourier
@reexport using ApproxFunOrthogonalPolynomials
@reexport using ApproxFunSingularities

import ApproxFunBase: Fun, UnsetSpace, VFun, UnivariateSpace, SumSpace, Space,
                    HeavisideSpace, PointSpace, IntervalOrSegment, ArraySpace,
                    TransformPlan, ITransformPlan, Evaluation,
                    Conversion, default_Fun, Derivative, Integral,
                    Dirichlet, domain, plan_transform,
                    plan_itransform, transform, domainspace,
                    rangespace, boundary, points, differentiate, integrate,
                    Segment, arclength, fromcanonical, checkpoints, evaluate,
                    coefficients, coefficientmatrix, clenshaw, ClenshawPlan,
                    SpaceOperator, InterlaceOperator, cfstype, pad!,
                    isapproxinteger, components, promotedomainspace, choosedomainspace,
                    AbstractProductSpace, MultivariateFun, BivariateSpace,
                    @calculus_operator, slnorm, sample, chop!, 𝒟, ∫, ⨜, ⨍

export ∫, ⨜, ⨍, 𝒟

import ApproxFunOrthogonalPolynomials: order

import BandedMatrices: bandwidths

import AbstractFFTs: Plan, fft, ifft
import FFTW: plan_fft, plan_ifft, plan_ifft!

import Base: convert, getindex, *, +, -, /, ^, \, sum, cumsum,
                first, last, isempty, zeros, promote_rule, real,
                # the following functions names are listed in Calculus.symbolic_derivatives_1arg(),
                # and methods are added to them here
                sqrt, cbrt, abs2, inv, log, log10, log2, log1p,
                exp, exp2, expm1, sin, cos, tan, sec, csc, cot,
                sind, cosd, tand, secd, cscd, cotd, asin, acos, atan, asec, acsc, acot,
                asind, acosd, atand, asecd, acscd, acotd, sinh, cosh, tanh, sech, csch,
                coth, asinh, acosh, atanh, asech, acsch, acoth, deg2rad, rad2deg

import LinearAlgebra: eigvals, dot, adjoint

import SpecialFunctions: erf, erfinv, erfc, erfcinv, erfi, gamma, lgamma, digamma, invdigamma,
                trigamma, airyai, airybi, airyaiprime, airybiprime, besselj0, besselj1,
                bessely0, bessely1, erfcx, dawson

# import Arpack: eigs


import FastTransforms: ChebyshevTransformPlan, plan_chebyshevtransform,
                        plan_chebyshevtransform!, plan_ichebyshevtransform,
                        plan_ichebyshevtransform!

using StaticArrays: SVector

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

##Testing
export bisectioninv


## Further extra features

include("Extras/Extras.jl")
include("Plot/Plot.jl")

#include("precompile.jl")
#_precompile_()

end #module
