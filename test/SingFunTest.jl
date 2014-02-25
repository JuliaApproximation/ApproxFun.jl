using ApproxFun


x=Fun(identity);
@assert sqrt(cos(π/2*x))[.1]-sqrt(cos(.1π/2)) < 10eps()