# Operators


```@meta
DocTestSetup = quote
    using ApproxFun
end
```

Linear operators between two spaces in ApproxFun are represented by
subtypes of `Operator`.  Every operator has a `domainspace` and `rangespace`.
That is, if a `Fun` `f` has the space `domainspace(op)`, then`op*f` is a
`Fun` with space `rangespace(op)`.

Note that the size of an operator is specified by the dimension of the domain
and range space.  

## Calculus operators

Differential and integral operators are perhaps the most useful type of operators
in mathematics.  Consider the derivative operator on `CosSpace`:
```jldoctest
julia> D = Derivative(CosSpace())
ConcreteDerivative:CosSpace(【0.0,6.283185307179586❫)→SinSpace(【0.0,6.283185307179586❫)
 0.0  -1.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    0.0  -2.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅     ⋅    0.0  -3.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅     ⋅     ⋅    0.0  -4.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅    0.0  -5.0    ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -6.0    ⋅     ⋅     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -7.0    ⋅     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -8.0    ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -9.0  ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  ⋱
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱

julia> f = Fun(θ->cos(cos(θ)), CosSpace());

julia> fp = D*f;

julia> fp(0.1) ≈ f'(0.1) ≈ sin(cos(0.1))*sin(0.1)
true
```

Here, we specified the domain space for the derivative operator, and it automatically
determined the range space:

```@meta
DocTestSetup = quote
    using ApproxFun
    D = Derivative(CosSpace())
    f = Fun(θ->cos(cos(θ)),CosSpace())
    fp = D*f
end
```

```jldoctest
julia> rangespace(D) == space(fp) == SinSpace()
true
```

Operators can be identified with infinite-dimensional matrices, whose entries
are given by the canonical bases in the domain and range space.  In this case,
the relevant formula is
$$D \cos k \theta = -k \sin k \theta.$$
That is, the `(k,k+1)`th entry is as follows:
```jldoctest
julia> k,j = 5,6;

julia> ej = Fun(domainspace(D),[zeros(j-1);1]);

julia> D[k,j] ≈ (D*ej).coefficients[k] ≈ -k
true
```

The `Chebyshev` space has the property that its derivatives are given by ultraspherical
spaces:
```jldoctest
julia> Derivative(Chebyshev())
ConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
```

## Functionals

A particularly useful class of operators are _functionals_, which map
from functions to scalar numbers.  These are represented by operators
of size `1 × ∞`: that is, infinite-dimensional analogues of row vectors.

As an example, the evaluation functional `f(0)`
on `CosSpace` has the form:

```jldoctest
julia> B = Evaluation(CosSpace(),0)
ConcreteEvaluation:CosSpace(【0.0,6.283185307179586❫)→ConstantSpace(Point(0))
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  ⋯

julia> B*f ≈ f(0)
true
```
As can be seen from the output, `rangespace(B)` is a `ConstantSpace(Point(0))`, a one-dimensional space used to represent scalars whose domain is a single point, `0`.


Closely related to functionals are operators with finite-dimensional range. For example,
the `Dirichlet` operator represents the restriction of a space to its boundary.
In the case, of `Chebyshev()`, this amounts to evaluation at the endpoints `±1`:
```jldoctest
julia> B = Dirichlet(Chebyshev())
ConcreteDirichlet:Chebyshev(【-1.0,1.0】)→2-element ArraySpace:
 ConstantSpace(Point(-1.0))
 ConstantSpace(Point(1.0))
 1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  ⋯
 1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  ⋯

julia> size(B)
(2, ∞)

julia> B*Fun(exp)
Fun(2-element ArraySpace:
 ConstantSpace(Point(-1.0))
 ConstantSpace(Point(1.0)) ,[0.367879, 2.71828])

julia> B*Fun(exp) ≈ Fun([exp(-1),exp(1)])
true
```

## Multiplication

A `Multiplication` operator sends a `Fun` to a `Fun` in the corresponding space by multiplying a given function. The `Multiplication` operators are presented in matrix form in `ApproxFun`.

```jldoctest
julia> x = Fun();

julia> M = Multiplication(1 + 2x + x^2, Chebyshev())
ConcreteMultiplication:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)
 1.5  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
 2.0  1.75  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
 0.5  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
  ⋅   0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅    ⋅
  ⋅    ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅    ⋅
  ⋅    ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅    ⋅
  ⋅    ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅    ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋱     ⋱    ⋱

julia> (M * x).coefficients == ((1 + 2x + x^2) * x).coefficients == M[1:4,1:2] * x.coefficients
true
```

It is possible for domain space and range space to be different under `Mulitplication`.

```jldoctest
julia> c = Fun(θ -> cos(θ), CosSpace());

julia> Multiplication(c, SinSpace())
ConcreteMultiplication:SinSpace(【0.0,6.283185307179586❫)→SinSpace(【0.0,6.283185307179586❫)
 8.974302782386682e-17  0.5                    …   ⋅                     ⋅
 0.5                    8.974302782386682e-17      ⋅                     ⋅
  ⋅                     0.5                        ⋅                     ⋅
  ⋅                      ⋅                         ⋅                     ⋅
  ⋅                      ⋅                         ⋅                     ⋅
  ⋅                      ⋅                     …   ⋅                     ⋅
  ⋅                      ⋅                         ⋅                     ⋅
  ⋅                      ⋅                         ⋅                     ⋅
  ⋅                      ⋅                        0.5                    ⋅
  ⋅                      ⋅                        8.974302782386682e-17  ⋱
  ⋅                      ⋅                     …   ⋱                     ⋱
```

If a function is given by the expansion
$$ f(\theta) = \sum_{n=1}^{\infty}  {f}_{n} * sin(n\theta) $$

Then the matrix above can be easily derived from
$$ cos(\theta) * f(\theta) = cos(\theta) \cdot (\sum_{n=1}^{\infty}  {f}_{n} \cdot sin(n\theta) $$
$$ = \sum_{n=1}^{\infty}  {f}_{n} \cdot cos(\theta) \cdot sin(n\theta) $$
$$ = \sum_{n=1}^{\infty}  {f}_{n} \cdot 0.5 \cdot ((sin(n-1)\theta) + (sin(n+1)\theta) $$
$$ = \sum_{n=1}^{\infty}  0.5 \cdot ({f}_{n-1} + {f}_{n+1}) \cdot sin(n\theta) $$.

## Algebraic manipulation of operators

Operators can be algebraically manipulated, provided that the domain and
range spaces are compatible, or can be made compatible.  
As a simple example, we can add the second derivative of a Fourier space to the
identity operator:
```jldoctest
julia> D2 = Derivative(Fourier(),2)
DerivativeWrapper:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)
 0.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅   -1.0    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅   -4.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅   -4.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅   -9.0    ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -9.0     ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -16.0     ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -16.0     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -25.0  ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱

julia> D2 + I
PlusOperator:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)
 1.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅   0.0   ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅   0.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅   -3.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅   -3.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅   -8.0    ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   -8.0     ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   -15.0     ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -15.0     ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -24.0  ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱
```

When the domain and range space are not the same, the identity operator
becomes a conversion operator.  That is, to represent `D+I` acting on the Chebyshev
space, we would do the following:

```jldoctest
julia> D = Derivative(Chebyshev())
ConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱

julia> C = Conversion(Chebyshev(),Ultraspherical(1))
ConcreteConversion:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 1.0  0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅   0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱


julia> D + C
PlusOperator:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 1.0  1.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅   0.5   2.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅    0.5   3.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅    0.5   4.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅    0.5   5.0  -0.5    ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   6.0  -0.5    ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   7.0  -0.5    ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   8.0  -0.5  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   9.0  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
```

ApproxFun can automatically determine the spaces, so if one writes
`D + I` it will translate it to `D + C`.  


Now consider the Fredholm integral operator of the second kind:

$$L u = u + {\rm e}^x \int_{-1}^1 u(x) {\rm d}x$$

We can construct this using

```jldoctest
julia> x = Fun();

julia> Σ = DefiniteIntegral(Chebyshev())
ConcreteDefiniteIntegral:Chebyshev(【-1.0,1.0】)→ConstantSpace
 2.0  0.0  -0.6666666666666666  0.0  …  0.0  -0.031746031746031744  0.0  ⋯

julia> L = I + exp(x)*Σ
LowRankPertOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)
 3.5321317555040164     0.0  …  -0.0401925675476828      0.0  ⋯
 2.260636415969941      1.0     -0.03588311771380859     0.0  ⋱
 0.5429906790681531     0.0     -0.008618899667748462    0.0  ⋱
 0.08867369969732766    0.0     -0.0014075190428147247   0.0  ⋱
 0.010948480884187475   0.0     -0.00017378541086011864  0.0  ⋱
 0.001085852623827888   0.0  …  -1.7235755933775998e-5   0.0  ⋱
 8.995464590859038e-5   0.0     -1.4278515223585775e-6   0.0  ⋱
 6.396872924803984e-6   0.0     -1.015376654730791e-7    0.0  ⋱
 3.9842496133455937e-7  0.0      0.9999999936757943      0.0  ⋱
 2.20735434510347e-8    0.0     -3.503737055719794e-10   1.0  ⋱
  ⋮                      ⋱   …    ⋱                       ⋱   ⋱

julia> u = cos(10x^2);

julia> (L*u)(0.1)
1.3777980523127336

julia> u(0.1) + exp(0.1)*sum(u)
1.3777980523127336
```

Note that `DefiniteIntegral` is a functional, i.e., a 1 × ∞ operator.  when
multiplied on the left by a function, it automatically constructs the operator
${\rm e}^x \int_{-1}^1 f(x) dx$
via

```@meta
DocTestSetup = quote
    using ApproxFun
    x = Fun()
    Q = DefiniteIntegral(Chebyshev())
end
```

```jldoctest
julia> M = Multiplication(exp(x),ConstantSpace())
ConcreteMultiplication:ConstantSpace→Chebyshev(【-1.0,1.0】)
 1.26607    
 1.13032    
 0.271495   
 0.0443368  
 0.00547424
 0.000542926
 4.49773e-5
 3.19844e-6
 1.99212e-7
 1.10368e-8
  ⋮         

julia> M*Σ
TimesOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)
 2.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯
 2.26064     0.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱
 0.542991    0.0  -0.180997     0.0     0.0  -0.0086189    0.0  ⋱
 0.0886737   0.0  -0.0295579    0.0     0.0  -0.00140752   0.0  ⋱
 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱
 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱
 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱
 6.39687e-6  0.0  -2.13229e-6   0.0     0.0  -1.01538e-7   0.0  ⋱
 3.98425e-7  0.0  -1.32808e-7   0.0     0.0  -6.32421e-9   0.0  ⋱
 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  0.0  ⋱
  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱
```
Note that `Q*exp(x)` applies the operator to a function.  To construct the operator that first multiplies by `exp(x)`, use `Q[exp(x)]`.  This is equivalent to `Q*Multiplication(exp(x),Chebyshev())`.


## Operators and space promotion

It is often more convenient to not specify a space explicitly, but rather
infer it when the operator is used.  For example, we can construct `Derivative()`,
which has the alias `𝒟`, and represents the first derivative on any space:
```jldoctest
julia> f = Fun(cos,Chebyshev(0..1)); (𝒟*f)(0.1)
-0.09983341664681705

julia> f = Fun(cos,Fourier()); (𝒟*f)(0.1)
-0.09983341664682804
```
Behind the scenes, `Derivative()` is equivalent to `Derivative(UnsetSpace(),1)`.
When multiplying a function `f`, the domain space is promoted before multiplying,
that is, `Derivative()*f` is equivalent to `Derivative(space(f))*f`.  

This promotion of the domain space happens even when operators have spaces attached.
This facilitates the following construction:
```jldoctest
julia> D = Derivative(Chebyshev());

julia> D^2
ConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(2,【-1.0,1.0】)
 ⋅  ⋅  4.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅   6.0   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅   8.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅   10.0    ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅   12.0    ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅   14.0    ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   16.0    ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   18.0  ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
```
Note that `rangespace(D) ≠ Chebyshev()`, hence the operators are not compatible.
Therefore, it has thrown away its domain space, and thus this is equivalent to
`Derivative(rangespace(D))*D`.





```@meta
DocTestSetup = nothing
```


## Concatenating operators

The concatenation functions `vcat`, `hcat` and `hvcat` are overriden for
operators to represent the resulting combined operator, now with
a `rangespace` or `domainspace` that is an `ArraySpace`.
