
# ApproxFun.jl NEWS

#### notes on release changes, ongoing development, and future planned work

## 0.5 (current master/development)

#### 0.5.0
- Drops support for Julia v0.5
- Uses IntervalSets.jl to support a..b
- Uses Padua points for `Chebyshev()^2` transform


## 0.4


#### 0.4.1
- `linsolve(A,b;kwds...)` -> `\(A,b;kwds...)`
- `transform(sp::Space,v,plan)` -> `plan*v`
- `PeriodicInterval()` now defaults to `PeriodicInterval(0,2π)`
- `points(::Chebyshev,n)` has reversed the order
- `Fun(cfs::Vector,sp::Space)` --> `Fun(sp::Space,cfs::Vector)`
- `Interval(a,b)` --> `Segment(a,b)` when `a` and `b` are not real-valued
- `Fun(f,[a,b])` --> `Fun(f,a..b)`, provided `a < b` are real-valued

#### 0.4.0
- Revamped PDE solving to use `qrfact`


