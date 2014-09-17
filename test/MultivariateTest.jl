using ApproxFun, Base.Test

u0   = TensorFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>Fun2D)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)[.1,.2]-sin(u0[.1,.2])|>abs < 10e-6


## Periodic
f=Fun2D((x,y)->cos(x)*sin(y),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(.1)*sin(.2)

f=Fun2D((x,y)->cos(cos(x)+sin(y)),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=points(f,1),y=points(f,2)]-values(f))<1000eps()

f=TensorFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval()^2)
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=points(f,1),y=points(f,2)]-values(f))<1000eps()

d=PeriodicInterval()^2
f=TensorFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
@test (f.'-f|>coefficients|>norm)< 10eps()



d=PeriodicInterval()^2
f=TensorFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 10000eps()

@test_approx_eq real(f)[.1,.2] f[.1,.2]
