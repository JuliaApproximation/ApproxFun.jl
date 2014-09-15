Pkg.add("PyPlot")
Pkg.add("RandomMatrices")

using ApproxFun

#The following samples eigenvalues of an n = 2 Gaussian Unitary Ensemble:

ff=(x,y)->(x.-y).^2.*exp(-x.^2/2.-y.^2/2);
f=Fun2D(ff,[-4.,4.],[-4.,4.]);
r=sample(f,5000)


#We can compare the histogram to the 1-point correlation

plot(sum(f,1)/sum(f)),PyPlot.plt.hist(vcat(r...),-4.:.1:4.,normed=true)
		
	
#We can compare the histograms of x with the GUE:

using RandomMatrices

plot(sum(f,1)/sum(f)),PyPlot.plt.hist(vcat([sqrt(2)eigvalrand(GaussianHermite(2),2) for k=1:10000]...),-4:.1:4,normed=true)