Funs is a package for approximating functions.  It currently supports intervals.  It is heavily influenced by the Matlab package chebfun (http://www.chebfun.org) and the Mathematica package RHPackage (http://www.maths.usyd.edu.au/u/olver/projects/RHPackage.html).

To construct an approximation of exp(x) on [-1,1], simply call

	f = Fun(exp,[-1,1])
	
The convention is to treat f as a vector.  In other words, we view functions as vectors with continuous indices.  To evaluate the function at a point, we use the vector notation:

	f[.1]
	
which will return a high accuracy approximation to  exp(.1).  Differentiation uses the diff operator:

	fp = diff(f)
	
Because the exponential is it's own derivative, the norm is small:

	norm(f-fp)
	
Similarly, cumsum defines an indefinite integration operator (not with a standard integration constant though!):

	intf = cumsum(f)
	intf = intf - intf[0] + f[0]
	norm(f - intf)
	
You can also add and multiply (use .* as these are treated like vectors):

	Fun(exp).*Fun(cos)				#gives a representation of exp(x)cos(x) on [-1,1]
	Fun(exp) + Fun(cos)				#gives a representation of exp(x) + cos(x) on [-1,1]	
	
Other operations including random number sampling using [Olver & Townsend 2013].  The following code samples 10,000 standard normals:

	f = Fun(x->exp(-x.^2),[-10,10])
	x = sample(f,10000)
	
We can apply this to any positive function.  The following samples the spectral density of an n = 4 Gaussian Unitary Ensemble):

	f = Fun(x->(9+72x.^2-192x.^4+512x.^6).*exp(-4x.^2),[-4,4])
	x = sample(f,40000)
		
		
Plotting is accomplished via Winston:

	using Winston
	plot(f)
	
We can compare the histograms of f with the GUE in RandomMatrices:



function sampleeigs(n,m)
    x=zeros(n,m);
    
    for k=1:m
        x[:,k]=eig(RandomMatrices.GaussianHermiteTridiagonalMatrix(n,2))[1];
    end
    
    reshape(x,n*m)
end

x2 = sampleeigs(4,10000)



	
	
References:
	
S. Olver & A. Townsend,  (2013), Fast inverse transform sampling in one and two dimensions, arXiv:1307.1223.