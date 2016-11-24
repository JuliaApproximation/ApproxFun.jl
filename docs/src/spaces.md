A `Space` is an abstract type that determines the basis (or frame) of an expansion.  Each `Domain` has a default space:

```julia
Space(Interval())     # returns Chebyshev(Interval())
```


##


The following commands must be overriden to provide features.


## Conversion of coefficients


`coefficients(cfs::Vector,fromspace::FunctionSpace,tospace::FunctionSpace)`

converts coefficients in `fromspace` to coefficients in `tospace`.  If not overriden, defaults to using `Conversion` operators (see below).

`canonicalspace(space::FunctionSpace)`

returns a space that is used as a default to implement missing functionality, e.g., evaluation.  Implement a Conversion operator or override `coefficients` to support this.



## Construction

`points(space::FunctionSpace,n::Integer)`

Return a grid.  Defaults to `points(domain(space),n)`

`transform(space::FunctionSpace,values::Vector)`

Transform values on the grid to coefficients. Defaults to `coefficients(transform(canonicalspace(space),values),canonicalspace(space),space)`


`itransform(space::FunctionSpace,coefficients::Vector)`

Transform coefficients back to values.  Defaults to using `canonicalspace` as in `transform`

## Evaluation

`evaluate(::FunctionSpace,coefficients::Vector,x::Number)`

Evaluates the expansion at a point x.  There are no requirements for the default behaviour of x not in the domain.  For example, Chebyshev will give the analytic continuation while Piecewise will return Nothing.  Defaults to converting to `canonicalspace` and evaluating.


## Helper routines

`spacescompatible`

Specifies equality of spaces while supporting `AnyDomain`

`Base.ones`

Return the `Fun` that represents the function one

## Conversion operator

`conversion_type(a::FunctionSpace,b::FunctionSpace)`

returns a `FunctionSpace` that has a banded conversion operator to both `a` and `b`.  Override `conversion_rule` when adding new `Conversion` operators.

`Conversion(a::FunctionSpace,b::FunctionSpace)`

while represent a conversion operator provided that `conversion_type` returns either `a` or `b`.  `addentries!` and `bandinds` must be overriden for this to function properly, see [[Operators|Operators]].  Otherwise, it will attempt to construct it via `Conversion(a,canonicalspace(a))*Conversion(canonicalspace(a),b)`
