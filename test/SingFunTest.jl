using ApproxFun


x=Fun(identity);
@assert sqrt(cos(π/2*x))[.1]-sqrt(cos(.1π/2)) < 10eps()


x=Fun(identity,[-2.,2.])
u=sqrt(4-x.^2)/(2π)
@assert abs(u[.1]-sqrt(4-.1^2)/(2π)) < 10eps()
@assert abs(sum(u) - 1) < 10eps()


f=Fun(x->x.*cot(π*x/2))
u=SingFun(f./(1-x.^2),1.,1.)
@assert abs(1./(.1.*cot(π*.1/2))-(1./u)[.1]) <10eps()

@assert abs((x./u)[.1]-tan(π*.1/2)) < 10eps()