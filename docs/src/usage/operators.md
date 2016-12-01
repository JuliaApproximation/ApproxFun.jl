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
 0.0  -1.0
       0.0  -2.0
             0.0  -3.0
                   0.0  -4.0
                         0.0  -5.0
                               0.0  -6.0
                                     0.0  -7.0
                                           0.0  -8.0
                                                 0.0  -9.0
                                                       0.0  ⋱
                                                            ⋱

julia> f = Fun(θ->cos(cos(θ)),CosSpace());

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
 0.0  1.0                                           
      0.0  2.0                                      
           0.0  3.0                                 
                0.0  4.0                            
                     0.0  5.0                       
                          0.0  6.0                  
                               0.0  7.0             
                                    0.0  8.0        
                                         0.0  9.0   
                                              0.0  ⋱
                                                   ⋱
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


## Algebraic manipulation of operators

Operators can be algebraically manipulated, provided that the domain and
range spaces are compatible, or can be made compatible.  
As a simple example, we can add the second derivative of a Fourier space to the
identity operator:
```jldoctest
julia> D2 = Derivative(Fourier(),2)
DerivativeWrapper:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)
 0.0   0.0                                                      
 0.0  -1.0   0.0                                                
       0.0  -1.0   0.0                                          
             0.0  -4.0   0.0                                    
                   0.0  -4.0   0.0                              
                         0.0  -9.0   0.0                        
                               0.0  -9.0    0.0                 
                                     0.0  -16.0    0.0          
                                            0.0  -16.0    0.0   
                                                   0.0  -25.0  ⋱
                                                           ⋱   ⋱

julia> D2 + I
PlusOperator:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)
 1.0  0.0                                                     
 0.0  0.0  0.0                                                
      0.0  0.0   0.0                                          
           0.0  -3.0   0.0                                    
                 0.0  -3.0   0.0                              
                       0.0  -8.0   0.0                        
                             0.0  -8.0    0.0                 
                                   0.0  -15.0    0.0          
                                          0.0  -15.0    0.0   
                                                 0.0  -24.0  ⋱
                                                         ⋱   ⋱
```

When the domain and range space are not the same, the identity operator
becomes a conversion operator.  That is, to represent `D+I` acting on the Chebyshev
space, we would do the following:

```jldoctest
julia> D = Derivative(Chebyshev())
ConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 0.0  1.0                                           
      0.0  2.0                                      
           0.0  3.0                                 
                0.0  4.0                            
                     0.0  5.0                       
                          0.0  6.0                  
                               0.0  7.0             
                                    0.0  8.0        
                                         0.0  9.0   
                                              0.0  ⋱
                                                   ⋱

julia> C = Conversion(Chebyshev(),Ultraspherical(1))
ConcreteConversion:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 1.0  0.0  -0.5                                             
      0.5   0.0  -0.5                                       
            0.5   0.0  -0.5                                 
                  0.5   0.0  -0.5                           
                        0.5   0.0  -0.5                     
                              0.5   0.0  -0.5               
                                    0.5   0.0  -0.5         
                                          0.5   0.0  -0.5   
                                                0.5   0.0  ⋱
                                                      0.5  ⋱
                                                           ⋱

julia> D + C
PlusOperator:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)
 1.0  1.0  -0.5                                             
      0.5   2.0  -0.5                                       
            0.5   3.0  -0.5                                 
                  0.5   4.0  -0.5                           
                        0.5   5.0  -0.5                     
                              0.5   6.0  -0.5               
                                    0.5   7.0  -0.5         
                                          0.5   8.0  -0.5   
                                                0.5   9.0  ⋱
                                                      0.5  ⋱
                                                           ⋱
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
 2.0  0.0  -0.666667  0.0  -0.133333  0.0  -0.0571429  0.0  -0.031746  0.0  ⋯

julia> L = I + exp(x)*Σ
LowRankPertOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)
 3.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯
 2.26064     1.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱
 0.542991    0.0   0.819003     0.0     0.0  -0.0086189    0.0  ⋱
 0.0886737   0.0  -0.0295579    1.0     0.0  -0.00140752   0.0  ⋱
 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱
 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱
 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱
 6.39687e-6  0.0  -2.13229e-6   0.0     1.0  -1.01538e-7   0.0  ⋱
 3.98425e-7  0.0  -1.32808e-7   0.0     0.0   1.0          0.0  ⋱
 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  1.0  ⋱
  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱

julia> u = cos(10x^2);

julia> (L*u)(0.1)
1.3777980523127333

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

julia> M*Q
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

## Operators and space promotion

It is often more convenient to not specify a space explicitly, but rather
infer it when the operator is used.  For example, we can construct `Derivative()`,
which has the alias `𝒟`, and represents the first derivative on any space:
```jldoctest
julia> f = Fun(cos,Chebyshev(0..1)); (𝒟*f)(0.1)
-0.09983341664681707

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
 0.0  0.0  4.0                                           
      0.0  0.0  6.0                                      
           0.0  0.0  8.0                                 
                0.0  0.0  10.0                           
                     0.0   0.0  12.0                     
                           0.0   0.0  14.0               
                                 0.0   0.0  16.0         
                                       0.0   0.0  18.0   
                                             0.0   0.0  ⋱
                                                   0.0  ⋱
                                                        ⋱
```
Note that `rangespace(D) ≠ Chebyshev()`, hence the operators are not compatible.
Therefore, it has thrown away its domain space, and thus this is equivalent to
`Derivative(rangespace(D))*D`.





```@meta
DocTestSetup = nothing
```
