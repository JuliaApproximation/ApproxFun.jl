using ApproxFun

#The following samples a 2D Cauchy Distribution
f = LowRankFun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),JacobiWeight(2.,2.,Line())^2)
@time x = ApproxFun.sample(f,100)
