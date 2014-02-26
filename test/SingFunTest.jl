using ApproxFun


x=Fun(identity);
@assert sqrt(cos(π/2*x))[.1]-sqrt(cos(.1π/2)) < 10eps()


x=Fun(identity,[-2.,2.])
u=sqrt(4-x.^2)/(2π)
@assert u[.1]-sqrt(4-.1^2)/(2π) < 10eps()
@assert abs(sum(u) - 1) < 10eps()
