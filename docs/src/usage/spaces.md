# Spaces

A `Space` is an abstract type whose subtypes indicate which space a function lives in.
This typically corresponds to the span of a (possibly infinite) basis.

## Classical orthogonal polynomial spaces

`Chebyshev`, `Ultraspherical`, `Jacobi`, `Hermite`, and `Laguerre` are spaces
corresponding to expansion in classical orthogonal polynomials.

Note that we always use the classical normalization: the basis are _not_ orthonormal.
This is because this normalization leads to rational recurrence relationships,
which are more efficient than their normalized counterparts. See
the [Digital Library of Mathematical Functions](https://dlmf.nist.gov/18) for
more information.


### Chebyshev space

The default space in ApproxFun is `Chebyshev`, which represents expansions
in Chebyshev polynomials:

$$f(x) = \sum_{k=0}^\infty f_k T_k(x)$$

where $T_k(x) = \cos k \,{\rm acos} x$, which are orthogonal polynomials
with respect to the weight
$$
{1 \over \sqrt{1-x^2}} \qquad\hbox{for}\qquad -1 \leq x \leq 1.
$$
Note that there is an intrinsic link between `Chebyshev` and `CosSpace`:  

$$g(\theta) = f(\cos \theta) = \sum_{k=0}^\infty f_k \cos k \theta$$

In other words:
```@meta
DocTestSetup = quote
    using ApproxFun
end
```
```jldoctest
julia> f=Fun(exp,Chebyshev());

julia> g=Fun(CosSpace(),f.coefficients); # specify the coefficients directly

julia> f(cos(0.1))
2.70473560723178

julia> g(0.1)
2.7047356072317794
```

### Ultraspherical spaces

A key tool for solving differential equations are the ultraspherical spaces,
encoded as `Ultraspherical(λ)` for `λ ≠ 0`,
which can be defined by the span of derivatives of Chebyshev polynomials, or alternatively
as polynomials orthogonal with respect to the weight
$$
(1-x^2)^{\lambda - 1/2} \qquad\hbox{for}\qquad -1 \leq x \leq 1.
$$

Note that `Ultraspherical(1)` corresponds to the Chebyshev basis of the
second kind: $U_k(x) = {\sin (k+1) {\rm acos} x \over \sin {\rm acos} x}$.  
The relationship with Chebyshev polynomials follows from trigonemetric
identities: $T_k'(x) = k U_{k-1}(x)$.  

Converting between ultraspherical polynomials (with integer orders) is
extremely efficient: it requires $O(n)$ operations, where $n$ is the number of coefficients.

### Jacobi spaces

`Jacobi(b,a)` represents the space spanned by the Jacobi polynomials, which
are orthogonal polynomials with respect to the weight
$$
(1+x)^b(1-x)^a
$$
using the standard normalization.


## Fourier and Laurent spaces

There are several different spaces to represent functions on periodic domains,
which are typically a `PeriodicInterval`, `Circle` or `PeriodicLine`.  

`CosSpace` represents expansion in cosine series:

$$f(\theta) = \sum_{k=0}^\infty f_k \cos k \theta$$

`SinSpace` represents expansion in sine series:

$$f(\theta) = \sum_{k=0}^\infty f_k \sin (k+1) \theta$$

`Taylor` represents expansion with only non-negative complex exponential terms:

$$f(\theta) = \sum_{k=0}^\infty f_k {\rm e}^{{\rm i} k \theta}$$

`Hardy{false}` represents expansion with only negative complex exponential terms:

$$f(\theta) = \sum_{k=0}^\infty f_k {\rm e}^{-{\rm i} (k+1) \theta}$$

`Fourier` represents functions that are sums of sines and cosines.  Note that
if a function has the form

$$f(\theta) = f_0 + \sum_{k=1}^\infty f_k^{\rm c} \cos k \theta + f_k^{\rm s} \sin k\theta$$

then the coefficients of the resulting `Fun` are order as $[f_0,f_1^{\rm s},f_1^{\rm c},…]$.
For example:

```jldoctest
julia> f = Fun(Fourier(),[1,2,3,4]);

julia> f(0.1)
4.979356652307978

julia> 1 + 2sin(0.1) + 3cos(0.1) + 4sin(2*0.1)
4.979356652307979
```

`Laurent` represents functions that are sums of complex exponentials.  Note that
if a function has the form

$$f(\theta) = \sum_{k=-\infty}^\infty f_k {\rm e}^{{\rm i} k \theta}$$

then the coefficients of the resulting `Fun` are order as $[f_0,f_{-1},f_1,…]$.
For example:

```jldoctest
julia> f = Fun(Laurent(),[1,2,3,4]);

julia> f(0.1)
9.895287137755096 - 0.694843906533417im

julia> 1 + 2exp(-im*0.1) + 3exp(im*0.1) + 4exp(-2im*0.1)
9.895287137755094 - 0.6948439065334167im
```


## Modifier spaces

Some spaces are built out of other spaces:

### `JacobiWeight`

 `JacobiWeight(β,α,space)`  weights `space`, which is typically `Chebyshev()` or `Jacobi(b,a)`,
 by a Jacobi weight `(1+x)^α*(1-x)^β`: in other words, if the basis for
`space` is $\psi_k(x)$ and the domain is the unit interval `-1 .. 1`, then the basis for
`JacobiWeight(β,α,space)` is $(1+x)^α(1-x)^β \psi_k(x)$. If the domain is
not the unit interval, then the basis is determined by mapping back to the
unit interval: that is, if $M(x)$ is the map dictated by `tocanonical(space, x)`,
where the canonical domain is the unit interval,
then the basis is $(1+M(x))^α(1-M(x))^β \psi_k(x)$. For example, if the domain is another
interval `a .. b`, then
$$
M(x) = {2x-b-a \over b-a}
$$
and the basis becomes
$$
\left({2 \over (b-a)}\right)^{\alpha+\eta}  (x-a)^α(b-x)^β \psi_k(x)
$$

### `SumSpace`

`SumSpace((space_1,space_2,…,space_n))` represents the direct sum of the spaces,
where evaluation is defined by adding up each component. A simple example is the
following, showing that the coefficients are stored by interlacing:
```jldoctest
julia> x = Fun(identity,-1..1);

julia> f = cos(x-0.1)*sqrt(1-x^2) + exp(x);

julia> space(f) # isa SumSpace
(1-x^2)^0.5[Chebyshev(【-1.0,1.0】)]⊕Chebyshev(【-1.0,1.0】)

julia> a, b = components(f);

julia> a(0.2) # cos(0.2-0.1)*sqrt(1-0.2^2)
0.9749009987500246

julia> b(0.2) # exp(0.2)
1.2214027581601699

julia> f(0.2) # a(0.2) + b(0.2)
2.1963037569101944

julia> norm(f.coefficients[1:2:end] - a.coefficients)
0.0

julia> norm(f.coefficients[2:2:end] - b.coefficients)
0.0
```
More complicated examples may interlace the coefficients using a different strategy.
Note that it is difficult to represent the first component of function $f$ by a Chebyshev series because the derivatives of $f$ at its boundaries blow up, whereas the derivative of a polynomial is a polynomial.

Note that `Fourier` and `Laurent` are currently implemented as `SumSpace`, but this
may change in the future.

### `PiecewiseSpace`

`PiecewiseSpace((space_1,space_2,…,space_n))` represents the direct sum of the spaces,
where evaluation is defined in a piecewise way. A simple example is the
following:
```jldoctest
julia> x = Fun(Domain(-1 .. 0) ∪ Domain( 1 .. 2));

julia> f = exp(x);

julia> a, b = components(f);

julia> f(-0.5) - a(-0.5)
0.0

julia> f(1.5) - b(1.5)
0.0

julia> f(0.5) # outside domain components
0.0

julia> norm(f.coefficients[2:2:end] - b.coefficients)
0.0

julia> norm(f.coefficients[1:2:end] - a.coefficients)
0.0
```
More complicated examples may interlace the coefficients using a different strategy.

### `ArraySpace`

`ArraySpace(::Array{<:Space})` represents the direct sum of the spaces,
where evaluation is defined in an array-wise manner. A simple example is the
following:
```jldoctest
julia> x = Fun(identity,-1..1);

julia> f = [exp(x); sqrt(1-x^2)*cos(x-0.1)];

julia> space(f)
2-element ArraySpace:
 Chebyshev(【-1.0,1.0】)             
 (1-x^2)^0.5[Chebyshev(【-1.0,1.0】)]

julia> a, b = components(f);

julia> norm(f(0.5) - [a(0.5); b(0.5)])
0.0

julia> norm(f.coefficients[1:2:end] - a.coefficients)
0.0

julia> norm(f.coefficients[2:2:end] - b.coefficients)
0.0
```
More complicated examples may interlace the coefficients using a different strategy.


### `TensorSpace`

`TensorSpace((space_1,space_2))` represents the tensor product of two spaces.
See documentation of `TensorSpace` for more details on how the coefficients
are interlaced. Note that more than two spaces is only partially supported.



## Unset space

`UnsetSpace` is a special space that is used as a stand in when a
space has not yet been determined, particularly by operators.  
