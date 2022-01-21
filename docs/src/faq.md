```@setup using-pkgs
using ApproxFun, Random
```

# Frequently Asked Questions

## Approximating functions

### How do I interpolate a function at a specified grid?

In the case where the grid is specified by `points(space,n)`, you can apply the default transform to data:

```@repl using-pkgs
S = Chebyshev(1..2);
p = points(S,20);  # the default grid
v = exp.(p);  # values at the default grid
f = Fun(S,ApproxFun.transform(S,v));
f(1.1)
exp(1.1)
```

ApproxFun has no inbuilt support for interpolating functions at other sets of points, but this can be accomplished manually by evaluating the basis at the set of points and using \:

```@repl using-pkgs
S = Chebyshev(1..2);
n = 50;
p = range(1,stop=2,length=n);  # a non-default grid
v = exp.(p);  # values at the non-default grid
V = Array{Float64}(undef,n,n);  # Create a Vandermonde matrix by evaluating the basis at the grid
for k = 1:n
   V[:,k] = Fun(S,[zeros(k-1);1]).(p)
end
f = Fun(S,V\v);
f(1.1)
exp(1.1)
```

Note that an evenly spaced grid suffers from instability for large `n`.  The easiest way around this is to use least squares with more points than coefficients, instead of interpolation:

```@repl using-pkgs
S = Chebyshev(1..2);
n = 100; m = 50;
p = range(1,stop=2,length=n);  # a non-default grid
v = exp.(p);  # values at the non-default grid
V = Array{Float64}(undef,n,m);  # Create a Vandermonde matrix by evaluating the basis at the grid
for k = 1:m
   V[:,k] = Fun(S,[zeros(k-1);1]).(p)
end
f = Fun(S,V\v);
f(1.1)
exp(1.1)
```

We can use this same approach for multivariate functions:

```@repl using-pkgs
S = Chebyshev(0..1)^2;
n = 1000; m = 50;
Random.seed!(0); x = rand(n); y = rand(n);
v = exp.(x .* cos.(y));  # values at the non-default grid
V = Array{Float64}(undef,n,m);  # Create a Vandermonde matrix by evaluating the basis at the grid
for k = 1:m
   V[:,k] = Fun(S,[zeros(k-1);1]).(x,y)
end
f = Fun(S,V\v);
f(0.1,0.2)
exp(0.1*cos(0.2))
```
