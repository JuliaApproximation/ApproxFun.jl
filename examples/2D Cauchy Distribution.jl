#The following samples a 2D Cauchy Distribution

	f = Fun2D((x,y)->1./(2π.*(x.^2 + y.^2 + 1).^(3/2),Line(),Line())
	x = sample(f,10000)
		
