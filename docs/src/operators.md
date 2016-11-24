
# Functionals

# Banded Operators

Add new BandedOperators

```julia
immutable MyBandedOperator{T} <: BandedOperator{T}
...
end
```

Override routines to work with built-in functionality like `*` and `\`.

```julia
domainspace(op::MyBandedOperator)
rangespace(op::MyBandedOperator)
```

Override to return the domain and range space of the operator:  The operator acts on coefficients in the domain space and returns coefficients in the range space.

```julia
bandinds(op::MyBandedOperator)
```

Override to return a tuple with the left and right bandwidths.

```julia
addentries!(op::MyBandedOperator,A,rows::Range)
```

Override to add the entries for the `rows` of `op` to the matrix data structure `A` (which may be a `BandedMatrix` or other `Matrix`-like data structure.)

# Calculus operators

Calculus operators `Derivative{S}` and `Integral{S}` are built-in concrete types that contains two fields: `domainspace::S` and `order`.  Overriding `bandinds(::Derivative{S})`, `rangespace(::Derivative{S})` and `addentries!(::Derivative{S},A,rows::Range)` to implement space-specific functionality.

(We omit extra template parameters that are not space dependent.)


# Multiplication operators

`Multiplication{S,V}` operators are a built-in concrete types that contains two fields: `f::Fun{S}` and `domainspace::V`.  Override `bandinds(::Multiplication{S,V})`, `addentries!(::Multiplication{S,V},A,rows::Range)` and `rangespace(::Multiplication{S,V})` to implement space-specific functionality.

# Conversion operators

`Conversion{S,V}` operators are a built-in concrete types that contains two fields: `domainspace::S` and `rangespace::V`.  Override `bandinds(:: Conversion{S,V})` and `addentries!(:: Conversion{S,V},A,rows::Range)` to implement space-specific functionality.


# Wrappers

`Derivative`, `Integral`, `Multiplication` and `Conversion` also have wrapper variants `DerivativeWrapper`, `IntegralWrapper`, `MultiplicationWrapper` and `ConversionWrapper` that allows for another operator to masquerade as a another operator.


# Space promotion

```julia
promotedomainspace(S::MyBandedOperator,domainspace::FunctionSpace)
promoterangespace(S::MyBandedOperator,rangespace::FunctionSpace)
```
can be overridden to change the domain/range space of the specified operator.  This allows lazily choosing the spaces by inferring them from other operators/the right-hand side.


```julia
choosedomainspace(S::MyBandedOperator,rangespace::FunctionSpace)
```
can be overridden to decide which domainspace should be chosen so that `MyBandedOperator` has the specified `rangespace`
