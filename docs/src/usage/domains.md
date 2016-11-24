`Domain` is an abstract type whose subtypes represent different domains.  Examples include `Interval`, `Ray`, `Line` and `Arc`.  Periodic domains include `PeriodicInterval`, `PeriodicLine` and `Circle`.  Vectors can be converted to domains:

```julia
Domain([-1,1])       # Returns Interval(-1.,1.)
Domain([0,Inf])      # Returns Ray(0.,true)
Domain([-Inf,Inf])   # Returns Line()
Domain([-Inf,-1,1])  # Returns Ray(-1.,false) ∪ Interval(0.,1.)
```

# Manipulating domains

Domains can be manipulated to make more complicated domains.  For example, you can take the union of an interval and a circle
```julia
Interval() ∪ Circle(3,0.5)    # equivalent to union(Interval(),Circle(3,0.5))
```
and the following creates a rectangle [0,1]^2:
```julia
rect=Interval(0,1)^2
```
Some other set operations are partially implemented:
```julia
Interval(0,2) ∩ Interval() # returns Interval(0,1)
```
