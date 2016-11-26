using ApproxFun, Base.Test
    import ApproxFun: testfunctional, testbandedbelowoperator, testbandedoperator

for S in (Fourier(Circle()),Laurent(Circle()),Taylor(Circle()),
            CosSpace(Circle()),JacobiWeight(-0.5,-0.5,Chebyshev()),
            JacobiWeight(-0.5,-0.5,Chebyshev(1.0..2.0+im)),
            JacobiWeight(0.5,0.5,Ultraspherical(1,1.0..2.0+im)))
    testfunctional(DefiniteLineIntegral(S))
end



# checks bug for
dom = Interval(-1.0,1) ∪ Interval(1.0+im,1+2im) ∪ Interval(-2.0+im,-1+2im)


⨍ = DefiniteLineIntegral(dom)
S = domainspace(⨍)
@test ApproxFun.bandinds(⨍) == (0,2)

f=Fun(S,rand(20))

@test_approx_eq DefiniteLineIntegral(dom[1])*f[1] linesum(f[1])
@test_approx_eq DefiniteLineIntegral(dom[2])*f[2] linesum(f[2])
@test_approx_eq DefiniteLineIntegral(dom[3])*f[3] linesum(f[3])

@test_approx_eq ⨍*f linesum(f)





#The first test checks the solution of the integral equation
# u(x) + \int_{-1}^{+1} \frac{e^{y} u(y)}{\sqrt{1-y^2}} dy = f
# on the interval -1..1.

x=Fun(identity)
w=1/sqrt(1-x^2)
d=domain(x)

Σ=DefiniteIntegral(d)
bandinds(Σ)

@test domainspace(Σ) ==
    JacobiWeight{Chebyshev{Interval{Float64}},Interval{Float64}}(-0.5,-0.5,Chebyshev())

L=I+Σ[exp(x)*w]
bandinds(L)
usol=sin(2x)
f=L*usol
u=L\f
@test norm(u-usol) <= 10eps()


#The second test checks the solution of the integro-differential equation
# u'(x) + x u(x) + \int_{-2}^{+2} sin(y-x) u(y) \sqrt{4-y^2} dy = f
# on the interval -2..2, with u(-2) = 1.

x=Fun(identity,-2..2)
w=sqrt(4-x^2)
d=domain(x)

D=Derivative(d)
B=ldirichlet(d)
Σ=DefiniteIntegral(.5,.5,d)

@test domainspace(Σ) ==
    JacobiWeight{Ultraspherical{Int,Interval{Float64}},Interval{Float64}}(.5,.5,Ultraspherical(1,d))

K=LowRankFun((x,y)->sin(y-x)*w(y),Ultraspherical(1,d),domainspace(Σ))

L=D+x+Σ[K]
usol=cospi(20x)
f=L*usol
u=[B;L]\[1.;f]


@test norm(u-usol) ≤ 200eps()


Σ = DefiniteIntegral()

f1=Fun(t->cos(cos(t)),-π..π)
f=Fun(t->cos(cos(t)),Laurent(-π..π))

@test_approx_eq sum(f1) Σ*f

f1=Fun(t->cos(cos(t))/t,Laurent(Circle()))
f2=Fun(t->cos(cos(t))/t,Fourier(Circle()))
@test_approx_eq Σ*f1 Σ*f2

f1=Fun(t->cos(cos(t)),Laurent(-π..π))
f2=Fun(t->cos(cos(t)),Fourier(-π..π))
@test_approx_eq Σ*f1 Σ*f2


## test over arcs


d=exp(im*Interval(0.1,0.2))
x=Fun(d)
w=1/(sqrt(abs(first(d)-x))*sqrt(abs(last(d)-x)))

@test_approx_eq linesum(w) DefiniteLineIntegral()*w


## Volterra integral equation

d = Interval(0.0,1.0)
V = Volterra(d)
K = LowRankFun((x,y)->sin(y-x),d^2)
L = I-V[K]

testbandedoperator(L)

f = Fun(exp,d)
@test domainspace(L) == Legendre(d)
@test rangespace(L) == Legendre(d)
@test bandrange(V) == -1:0
u = L\f
@test norm(L*u-f) ≤ 20eps()



## Check DefiniteIntegral

for S in (JacobiWeight(0.5,0.5,Ultraspherical(1,Interval(-2,-1))),
          JacobiWeight(0.5,0.5,Ultraspherical(1,Interval(-2,-1+2im))),
          JacobiWeight(1.5,1.5,Ultraspherical(2,Interval(-2,-1+2im))),
          JacobiWeight(-0.5,-0.5,Chebyshev(Interval(-2,-1+2im))),
          JacobiWeight(0.67,0.123,Chebyshev(Interval(-2,-1+2im))),
          JacobiWeight(0.67,0.123,Ultraspherical(1,Interval(-2,-1+2im))))
    f=Fun(S,[1.,2.,3.])
    @test_approx_eq DefiniteIntegral(space(f))*f sum(f)
    @test_approx_eq DefiniteLineIntegral(space(f))*f linesum(f)
end



## Fredholm integral

K=LowRankFun((x,y)->cos(x-y),Interval()^2)
Σ=DefiniteIntegral(Chebyshev())
testbandedbelowoperator(Σ[K])
