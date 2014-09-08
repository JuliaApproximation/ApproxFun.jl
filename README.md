`ApproxFun` is a package for approximating functions.  It currently supports intervals, the real 
line, periodic intervals and the unit circle.  It is heavily influenced by the Matlab 
package [`chebfun`](http://www.chebfun.org) and the Mathematica package [`RHPackage`](http://www.maths.usyd.edu.au/u/olver/projects/RHPackage.html).




To construct an approximation of `exp(x)` on `[-1,1]`, call


    using ApproxFun
	f = Fun(exp,[-1,1])
	
The convention is to treat `f` as a vector.  In other words, we view functions as vectors 
with continuous indices.  To evaluate the function at a point, we use the vector notation:

	f[.1]
	
which will return a high accuracy approximation to `exp(.1)`.  


# Differentiation and integration	


Differentiation uses the `diff` operator:

	fp = diff(f)
	
Because the exponential is its own derivative, the `norm` is small:

	norm(f-fp)
	
Similarly, `cumsum` defines an indefinite integration operator:

	g = cumsum(f)
	g = g + f[-1]
	norm(f - g)
	
You can also add and multiply (use `.*` as these are treated like vectors):

	Fun(exp).*Fun(cos)				#gives a representation of exp(x)cos(x) on [-1,1]
	Fun(exp) + Fun(cos)				#gives a representation of exp(x) + cos(x) on [-1,1]	
	
# Sampling	

Other operations including random number sampling using [Olver & Townsend 2013].  The 
following code samples 10,000 standard normals:

	f = Fun(x->exp(-x.^2),[-10,10])
	x = ApproxFun.sample(f,10000)
    plot(f)             # 2D plotting requires Gadfly
	Gadfly.plot(x=x,Gadfly.Geom.histogram)
	
We can apply this to any positive smooth PDF.  



There is also support for Fourier representations of functions on periodic intervals.  
Use `FFun` to ensure that the representation is periodic:

	f = FFun(cos)
	plot(f)

The default domain is `[-π,π]`.  



Differentiation is again accomplished with `diff`, so the following will be small:

	norm(diff(f) + FFun(sin))

Indefinite integration is only supported when the zeroth Fourier coefficient is zero:
	
	norm(cumsum(f) - FFun(sin))	

	
	
Alternatively, a Laurent series can be constructed on the unit circle:

	c = FFun(cos,Circle())
	


# Solving ordinary differential equations


We can solve ODEs, the following solves the Airy equation `u' = x u` as a BVP on `[-1000,10]`:


	x=Fun(identity,[-1000.,15.])
   	d=x.domain
	D=diff(d)
	u = [dirichlet(d),D^2 - x] \ [airyai(d.a),0.]
	
	plot(u)
	
# Solving partial differential equations

We can solve PDEs, the following solves Helmholtz `Δu + 100u=0` with `u(±1,y)=u(x,±1)=1`


    d=Interval()⊗Interval()
    
    u=[dirichlet(d),lap(d)+100I]\ones(4);
    contour(u)
	
	
# References

S. Olver & A. Townsend (2013), A fast and well-conditioned spectral method, SIAM Review, 55:462Ð489.
	
S. Olver & A. Townsend  (2013), Fast inverse transform sampling in one and two dimensions, arXiv:1307.1223.

A. Townsend & S. Olver (2014), The automatic solution of partial differential equations using a global spectral method

