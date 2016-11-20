# Frequently Asked Questions

## Approximating functions

### How do I interpolate a function at a specified grid?

In the case where the grid is specified by `points(space,n)`, you can apply the default transform to data:

```julia
S = Chebyshev([1,2])  
p = points(S,20) # the default grid
v = exp.(p)      # values at the default grid
f = Fun(ApproxFun.transform(S,vals),S)
```

ApproxFun has no inbuilt support for interpolating functions at other sets of points, but this can be accomplished manually by evaluating the basis at the set of points and using \:

```julia
S = Chebyshev([1,2])  
n = 50
p = linspace(1,2,n)  # a non-default grid
v = exp.(p)           # values at the non-default grid
# Create a Vandermonde matrix by evaluating the basis at the grid
V = Array(Float64,n,n)
for k = 1:n
    V[:,k] = Fun([zeros(k-1);1],S)(p)
end
f = Fun(V\v,S)   
```

Note that an evenly spaced grid suffers from instability for large `n`.  The easiest way around this is to use least squares with more points than coefficients, instead of interpolation:

```julia
S = Chebyshev([1,2])  
n = 100; m = 50
p = linspace(1,2,n)  # a non-default grid
v = exp.(p)           # values at the non-default grid
# Create a Vandermonde matrix by evaluating the basis at the grid
V = Array(Float64,n,m)
for k = 1:m
    V[:,k] = Fun([zeros(k-1);1],S)(p)
end
f = Fun(V\v,S)   
```
