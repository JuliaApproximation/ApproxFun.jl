# Domains

`Domain` is an abstract type whose subtypes represent oriented domains on which we wish to approximate functions.  Examples include `Interval`, `Ray`, `Line` and `Arc`.  Periodic domains include `PeriodicSegment`, `PeriodicLine` and `Circle`.

## Relationship with spaces

Every domain `d` has a default space, constructed via `Space(d)`.  For example, the default space for `ChebyshevInterval()` is `Chebyshev(ChebyshevInterval())`, which is efficient for representing smooth functions.  On the other hand, the default space for `PeriodicSegment()` is `Fourier(PeriodicSegment())`, which uses trigonometric polynomials to approximate periodic functions.  

## Manipulating domains

Domains can be manipulated to make more complicated domains.  For example, you can take the union of an interval and a circle

```julia
ChebyshevInterval() ∪ Circle(3,0.5)  # equivalent to union(ChebyshevInterval(),Circle(3,0.5))
```
and the following creates a rectangle `[0,1]^2`:

```julia
rect=Interval(0,1)^2
```

Some other set operations are partially implemented:

```julia
Interval(0,2) ∩ ChebyshevInterval()  # returns Interval(0,1)
```
