using ApproxFun

x=Fun()
u0=0.0x
N=u->[u(-1.),u(1.),0.001u''+6*(1-x^2)*u'+u^2-1.]
u=newton(N,u0)


x=Fun()
u0=0.0x
N=u->[u(-1.),u(1.),0.1u''+sin(u)-1.]
u=newton(N,u0)
