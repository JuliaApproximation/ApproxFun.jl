using ApproxFun

x=Fun()
d=domain(x)

N=u->[u[-1.],u[1.],u''+2*(1-x^2)*u'+u^2-1.]


@time u=newton(N,zeros(d))
ApproxFun.plot(u)

