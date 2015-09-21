using ApproxFun, Base.Test






## Try constructor variants

ff=(x,y)->exp(-10(x+.2)^2-20(y-.1)^2)*cos(x*y)
gg=x->exp(-10(x[1]+.2)^2-20(x[1]-.1)^2)*cos(x[1]*x[2])
f=Fun(ff,Interval()^2,10000)
@test_approx_eq f(0.,0.) ff(0.,0.)

f=Fun(gg,Interval()^2,10000)
@test_approx_eq f(0.,0.) ff(0.,0.)

f=Fun(ff,Interval()^2)
@test_approx_eq f(0.,0.) ff(0.,0.)
f=Fun(gg,Interval()^2)
@test_approx_eq f(0.,0.) ff(0.,0.)


f=Fun(ff)
@test_approx_eq f(0.,0.) ff(0.,0.)
f=Fun(gg)
@test_approx_eq f(0.,0.) ff(0.,0.)



## ProductFun
u0   = ProductFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>LowRankFun)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)(.1,.2)-sin(u0(.1,.2))|>abs < 10e-4

## LowRankFun

F = LowRankFun((x,y)->besselj0(10(y-x)),Chebyshev(),Chebyshev())

@test_approx_eq F(.123,.456) besselj0(10(.456-.123))

G = LowRankFun((x,y)->besselj0(10(y-x));method=:Cholesky)

@test_approx_eq G(.357,.246) besselj0(10(.246-.357))

F = LowRankFun((x,y)->hankelh1(0,10abs(y-x)),Chebyshev([1.0,2.0]),Chebyshev([1.0im,2.0im]))

@test_approx_eq F(1.5,1.5im) hankelh1(0,10abs(1.5im-1.5))




## Periodic
f=LowRankFun((x,y)->cos(x)*sin(y),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f(.1,.2) cos(.1)*sin(.2)

f=LowRankFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f(.1,.2) cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=ApproxFun.vecpoints(f,1),y=ApproxFun.vecpoints(f,2)]-values(f))<10000eps()

f=ProductFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval()^2)
@test_approx_eq f(.1,.2) cos(cos(.1)+sin(.2))
x,y=points(f)
@test norm(Float64[cos(cos(x[k,j])+sin(y[k,j])) for k=1:size(f,1),j=1:size(f,2)]-values(f))<10000eps()

d=PeriodicInterval()^2
f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
@test (f.'-f|>coefficients|>norm)< 1000eps()







## Functional*Fun

d=Interval()
B=ldirichlet(d)
f=ProductFun((x,y)->cos(cos(x)*sin(y)),d^2)
@test norm(B*f-Fun(y->cos(cos(-1)*sin(y)),d))<20000eps()
@test norm(f*B-Fun(x->cos(cos(x)*sin(-1)),d))<20000eps()
