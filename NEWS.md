
# ApproxFun.jl NEWS

### notes on release changes, ongoing development, and future planned work

#### 0.11
 - Split into ApproxFunBase.jl, ApproxFunFourier.jl, ApproxFunOrthogonalPolynomials.jl, and ApproxFunSingularities.jl

#### 0.10.4
 - Update for latest BandedMatrices.jl
 - Support for symmetric discretizations
 - Bug fixes

#### 0.10.3
 - Bug fixes
 - Update for latest BlockBandedMatrices.jl

#### 0.10.2
 - Bug fixes for special functions on large intervals

#### 0.10.1 
 - Update for latest BlockBandedMatrices.jl
 
#### 0.10.0
 - Use DomainSets.jl for representing domains
 - Use InfiniteArrays.jl for representing block sizes
 - Improve support for Laguerre
 - Remove `bandinds` in favour of `bandwidths`
 - Fixes for Hermite polynomials
 - Support latest BandedMatrices and BlockBandedMatrices
  
#### 0.9.0
 - Upgrade to Julia v1.0
 - Improve support for Heaviside functions

#### 0.8.2
 - Support integrating `WeightedLaguerre`
 - Prepare for Julia v0.7

#### 0.8.1 
 - Support BandedMatrices v0.5
 - Bug fixes and minor improvements

#### 0.8.0
- In-place transform functionality for `Fourier`
- Use BandedBlockBandedMatrices.jl for more reliable solution of PDEs
- Adds `jumplocations` for finding discontinuities of a piecewise `Fun` (thanks lcw)
- Adds `KroneckerDelta` (thaks marcusdavidwebb)
- Improvements for `BigFloat` with `JacobiWeight`
- Fix several bugs and performance enhancements

#### 0.7.1
- `F` was renamed `DFunction` for dynamic funtion
- `Fun`s are now subtypes of `Function`
- Support `f^k` for polynomial `f` and real `k`
- Fix several bugs and performance enhancements

#### 0.7.0
- Examples moved to [ApproxFunExamples](https://github.com/JuliaApproximation/ApproxFunExamples) repository
- `vcat`, `hcat` and `hvcat` of operators now returns an `Operator`
- `Dirichlet` and `Neumann` operators supported in 1D, replacing
`dirichlet` and `neumann`
- Improved broadcasting over a `Fun`
- Add multivariate `DefiniteIntegral`
- `Evaluation(first)`/`Evaluation(last)` now used for
evaluating at the boundary of an interval
- Fast evaluation for `SumSpace` and `PiecewiseSpace`
- Support general `AbstractVector` coefficients, to allow for sub-views
- Simplifies `Space` and `Domain` templated variables
- No longer abuses vector notation for `PiecewiseSpace`, `SumSpace` and `TensorSpace`
- Drops support for Julia v0.5

#### 0.6.1
- Uses more sophisticated chopping algorithm from [Aurentz & Trefethen 2016](https://people.maths.ox.ac.uk/trefethen/aurentz_trefethen_toms_final.pdf)

#### 0.6.0
- Adds support for Julia v0.6
- Replaces FixedSizeArrays.jl dependancy with StaticArrays.jl
- Auto-vectorization `f([1,2,3])` is removed in favour of broadcasting `f.([1,2,3])`


#### 0.5.0
- Drops support for Julia v0.4
- Uses IntervalSets.jl to support a..b
- Uses Padua points for `Chebyshev()^2` transform


#### 0.4.1
- `linsolve(A,b;kwds...)` -> `\(A,b;kwds...)`
- `transform(sp::Space,v,plan)` -> `plan*v`
- `PeriodicSegment()` now defaults to `PeriodicSegment(0,2π)`
- `points(::Chebyshev,n)` has reversed the order
- `Fun(cfs::Vector,sp::Space)` --> `Fun(sp::Space,cfs::Vector)`
- `Interval(a,b)` --> `Segment(a,b)` when `a` and `b` are not real-valued
- `Fun(f,[a,b])` --> `Fun(f,a..b)`, provided `a < b` are real-valued

#### 0.4.0
- Revamped PDE solving to use `qrfact`
