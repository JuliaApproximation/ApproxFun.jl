using ApproxFun

#The following samples a 2D Cauchy Distribution

f = Fun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),Line()^2)
x = ApproxFun.sample(f,100)
