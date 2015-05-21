[![Build Status](https://travis-ci.org/ApproxFun/ApproxFun.jl.svg?branch=master)](https://travis-ci.org/ApproxFun/ApproxFun.jl) [![Coverage Status](https://img.shields.io/coveralls/ApproxFun/ApproxFun.jl.svg)](https://coveralls.io/r/ApproxFun/ApproxFun.jl?branch=master)

`ApproxFun` is a package for approximating functions. It is heavily influenced by the Matlab 
package [`Chebfun`](http://www.chebfun.org) and the Mathematica package [`RHPackage`](http://www.maths.usyd.edu.au/u/olver/projects/RHPackage.html).



Take your two favourite functions on an interval and create approximations to them as simply as:

```julia
using ApproxFun
x = Fun(identity,[0.,10.])
f = sin(x^2)
g = cos(x)
```

To evaluate functions at a point, we use the vector notation, and `f[.1]` will return a high
accuracy approximation to `sin(0.01)`. All the algebraic manipulations of functions 
are supported and more.  For example, we can add `f` and `g^2` together and compute 
the roots and extrema:

```julia
h = f + g^2
r = roots(h)
rp = roots(differentiate(h))
ApproxFun.plot(h)                      # using PyPlot
PyPlot.plot(r,h[r],"og",rp,h[rp],"or") # using PyPlot
```

![Extrema](https://github.com/ApproxFun/ApproxFun.jl/raw/master/images/extrema.png)


# Differentiation and integration	


Notice from above that to find the extrema, we used the `differentiate` operator. Several other `Julia`
base functions are overloaded for the purposes of calculus. Because the exponential is its own
derivative, the `norm` is small:

```julia
f = Fun(x->exp(x),[-1.,1.])
fp = differentiate(f)
norm(f-fp)
```

Similarly, `cumsum` defines an indefinite integration operator:

```julia
g = cumsum(f)
g = g + f[-1]
norm(f-g)
```

`Fun`s in `ApproxFun` are instances of `Julia` types with one field to store coefficients and another
to describe the function space. Similarly, each function space has one field describing 
its domain. Let's explore:

```julia
x = Fun(identity)
f = exp(x)
g = f/sqrt(1-x^2)
space(f)
space(g)
```

In this case, `f` is in the `Ultraspherical{0}` space on the domain `Interval(-1.0,1.0)`, and
`g` is in the decorated `JacobiWeight{Ultraspherical{0}}` space. The absolute value is 
another case where space promotion is inferred from the operation:

```julia
f = Fun(x->cospi(5x))
g = abs(f)
space(f)
space(g)
```

Algebraic and differential operations are also implemented where possible.


# Solving ordinary differential equations


The following solves the Airy ODE `u'' - x u = 0` as a BVP on `[-1000,200]`:

```julia
x = Fun(identity,[-1000.,200.])
d = domain(x)
D = Derivative(d)
B = dirichlet(d)
L = D^2 - x
u = [B;L] \ [airyai(d.a);airyai(d.b)]
ApproxFun.plot(u)						    # Requires Gadfly or PyPlot
```

![Airy](https://github.com/ApproxFun/ApproxFun.jl/raw/master/images/airy.png)


# Periodic functions


There is also support for Fourier representations of functions on periodic intervals. 
Specify the space `Fourier` to ensure that the representation is periodic:

```julia
f = Fun(cos,Fourier([-π,π]))
norm(differentiate(f) + Fun(sin,Fourier([-π,π]))
```

Due to the periodicity, Fourier representations allow for the asymptotic savings of `2/π` 
in the number of coefficients that need to be stored compared with a Chebyshev representation. 
ODEs can also be solved when the solution is periodic:

```julia
s = Chebyshev([-π,π])
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
B = periodic(s,0)
uChebyshev = [B;L]\[0.;f]

s = Fourier([-π,π])
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
uFourier = L\f

length(uFourier)/length(uChebyshev),2/π
ApproxFun.plot(uFourier)						    # Requires Gadfly or PyPlot
```

![Periodic](https://github.com/ApproxFun/ApproxFun.jl/raw/master/images/periodic.png)


# Sampling	


Other operations including random number sampling using [Olver & Townsend 2013].  The 
following code samples 10,000 from a PDF given as the absolute value of the sine function on `[-5,5]`:

```julia
f = abs(Fun(sin,[-5,5]))
x = ApproxFun.sample(f,10000)
ApproxFun.plot(f/sum(f))                           # Requires Gadfly or PyPlot
PyPlot.plt.hist(x;normed=true,bins=[-5.:.1:5.])
```

![Sampling](https://github.com/ApproxFun/ApproxFun.jl/raw/master/images/sample.png)


# Solving partial differential equations


We can solve PDEs, the following solves Helmholtz `Δu + 100u=0` with `u(±1,y)=u(x,±1)=1`
on a square

```julia
d = Interval()^2          					# Defines a rectangle

u = [dirichlet(d);lap(d)+100I]\ones(4)		# First four entries of rhs are 
    											# boundary conditions
ApproxFun.contour(u)						# Requires Gadfly or PyPlot
```

The following solves Poisson `Δu =f` with zero Dirichlet conditions
on a disk

```julia
d = Disk()
f = Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),d) 
u = [dirichlet(d);lap(d)]\Any[0.,f]
ApproxFun.plot(u)                           # Requires Gadfly or PyPlot
```

We can also evolve PDEs.  The following solves advection—diffusion 
`u_t = 0.01Δu - 4u_x -3u_y` on a rectangle

```julia
d = Interval()^2
u0 = Fun((x,y)->exp(-40(x-.1)^2-40(y+.2)^2),d)
B = dirichlet(d)
D = Derivative(Interval())
L = (0.01D^2-4D)⊗I + I⊗(0.01D^2-3D)
h = 0.002
timeevolution(B,L,u0,h)                    # Requires GLPlot
```

The following solves beam equation `u_tt + Δ^2u = 0`
on a disk

```julia
d = Disk()
u0 = Fun((x,y)->exp(-50x^2-40(y-.1)^2)+.5exp(-30(x+.5)^2-40(y+.2)^2),d)
B= [dirichlet(d),neumann(d)]
L = -lap(d)^2
h = 0.001
timeevolution(2,B,L,u0,h)                 # Requires GLPlot
```



	
# References

S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, arXiv:1409.5529, to appear in HPTCDL 2014

A. Townsend & S. Olver (2014), The automatic solution of partial differential equations using a global spectral method, arXiv:1409:2789

S. Olver & A. Townsend (2013), Fast inverse transform sampling in one and two dimensions, arXiv:1307.1223

S. Olver & A. Townsend (2013), A fast and well-conditioned spectral method, SIAM Review, 55:462–489
	
