# Frequently Asked Questions

## Approximating functions

### How do I interpolate a function at a specified grid?

In the case where the grid is specified by `points(space,n)`, you can apply the default transform to data:

```@meta
DocTestSetup = quote
    using ApproxFun
end
```

```jldoctest
julia> S = Chebyshev(1..2);

julia> p = points(S,20); # the default grid

julia> v = exp.(p);      # values at the default grid

julia> f = Fun(S,ApproxFun.transform(S,v));

julia> f(1.1)
3.0041660239464347

julia> exp(1.1)
3.0041660239464334
```

ApproxFun has no inbuilt support for interpolating functions at other sets of points, but this can be accomplished manually by evaluating the basis at the set of points and using \:

```jldoctest
julia> S = Chebyshev(1..2);

julia> n = 50;

julia> p = linspace(1,2,n);   # a non-default grid

julia> v = exp.(p);           # values at the non-default grid

julia> V = Array(Float64,n,n); # Create a Vandermonde matrix by evaluating the basis at the grid

julia> for k = 1:n
           V[:,k] = Fun(S,[zeros(k-1);1]).(p)
       end

julia> f = Fun(S,V\v);

julia> f(1.1)
3.0041660228311926

julia> exp(1.1)
3.0041660239464334
```

Note that an evenly spaced grid suffers from instability for large `n`.  The easiest way around this is to use least squares with more points than coefficients, instead of interpolation:

```jldoctest
julia> S = Chebyshev(1..2);

julia> n = 100; m = 50;

julia> p = linspace(1,2,n);   # a non-default grid

julia> v = exp.(p);           # values at the non-default grid

julia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid

julia> for k = 1:m
           V[:,k] = Fun(S,[zeros(k-1);1]).(p)
       end

julia> f = Fun(S,V\v);

julia> f(1.1)
3.004166023946434

julia> exp(1.1)
3.0041660239464334
```

We can use this same approach for multivariate functions:

```jldoctest
julia> S = Chebyshev(0..1)^2;

julia> n = 1000; m = 50;

julia> srand(0); x = rand(n); y = rand(n);

julia> v = exp.(x .* cos(y));  # values at the non-default grid

julia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid

julia> for k = 1:m
          V[:,k] = Fun(S,[zeros(k-1);1]).(x,y)
       end


julia> f = Fun(S,V\v);

julia> f(0.1,0.2)
1.1029700685084018

julia> exp(0.1*cos(0.2))
1.1029701284210731
```


```@meta
DocTestSetup = nothing
```
