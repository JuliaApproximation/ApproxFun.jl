#The following samples the spectral density of an n = 4 Gaussian Unitary Ensemble):

	f = Fun(x->(9+72x.^2-192x.^4+512x.^6).*exp(-4x.^2),[-4,4])
	x = sample(f,40000)
		
		
# The PDF compared to the samples

    plot(f/sum(f));PyPlot.plt.hist(x,100,normed=true)
	
#We can compare the histograms of x with the GUE:

    using RandomMatrices

    PySide.plt.hist(vcat([eigvalrand(GaussianHermite(2),2) for k=1:1000]...),-4:.1:4)