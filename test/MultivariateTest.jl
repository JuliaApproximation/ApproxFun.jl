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


# Fun +-* constant
f=Fun((x,y)->exp(x)*cos(y))

@test_approx_eq f(0.1,0.2)+2 (f+2)(0.1,0.2)
@test_approx_eq f(0.1,0.2)-2 (f-2)(0.1,0.2)
@test_approx_eq f(0.1,0.2)*2 (f*2)(0.1,0.2)



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



## 1D in 2D

d=Interval((0.,0.),(1.,1.))
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))
@test_approx_eq f(FixedSizeArrays.Vec(0.5,0.5)) exp(-0.5-2cos(0.5))

f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d,20)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d,20)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))


d=Circle((0.,0.),1.)
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),Fourier(d),40)
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))
@test_approx_eq f(FixedSizeArrays.Vec(cos(0.1),sin(0.1))) exp(-cos(0.1)-2cos(sin(0.1)))

f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d),40)
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))


f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d))
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))

println("    Calculus tests")

## Sum

ff=(x,y)->(x-y)^2*exp(-x^2/2.-y^2/2)
f=Fun(ff,[-4.,4.],[-4.,4.])


@test_approx_eq sum(f,1)(0.1) 2.5162377980828357
f=LowRankFun(f)
@test_approx_eq evaluate(f.A,0.1) map(f->f(0.1),f.A)


# 2d derivative (issue #346)
let d = Chebyshev()^2
    f = Fun((x,y) -> sin(x) * cos(y), d)
    C=Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical{1}()⊗Ultraspherical{1}())
    @test_approx_eq (C*f)(0.1,0.2) f(0.1,0.2)
    Dx = Derivative(d, [1,0])
    f = Fun((x,y) -> sin(x) * cos(y), d)
    fx = Fun((x,y) -> cos(x) * cos(y), d)
    @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
    Dy = Derivative(d, [0,1])
    fy = Fun((x,y) -> -sin(x) * sin(y), d)
    @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
    L=Dx+Dy
    @test_approx_eq (L*f)(0.2,0.3) (fx(0.2,0.3)+fy(0.2,0.3))

    B=ldirichlet(d[1])⊗ldirichlet(d[2])
    @test_approx_eq B*f f(-1.,-1.)

    B=Evaluation(d[1],0.1)⊗ldirichlet(d[2])
    @test_approx_eq B*f f(0.1,-1.)

    B=Evaluation(d[1],0.1)⊗Evaluation(d[2],0.3)
    @test_approx_eq B*f f(0.1,0.3)

    B=Evaluation(d,(0.1,0.3))
    @test_approx_eq B*f f(0.1,0.3)
end

let d = Space([0,1]) * Space([0,2])
    Dx = Derivative(d, [1,0])
    f = Fun((x,y) -> sin(x) * cos(y), d)
    fx = Fun((x,y) -> cos(x) * cos(y), d)
    @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
    Dy = Derivative(d, [0,1])
    fy = Fun((x,y) -> -sin(x) * sin(y), d)
    @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
    L=Dx+Dy
    @test_approx_eq (L*f)(0.2,0.3) (fx(0.2,0.3)+fy(0.2,0.3))

    B=ldirichlet(d[1])⊗ldirichlet(d[2])
    @test abs(B*f-f(0.,0.)) ≤ 10eps()

    B=Evaluation(d[1],0.1)⊗ldirichlet(d[2])
    @test_approx_eq B*f f(0.1,0.)

    B=Evaluation(d[1],0.1)⊗Evaluation(d[2],0.3)
    @test_approx_eq B*f f(0.1,0.3)

    B=Evaluation(d,(0.1,0.3))
    @test_approx_eq B*f f(0.1,0.3)
end
