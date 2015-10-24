
## Intro

using ApproxFun,Base.Test
x = Fun(identity,[0.,10.])
f = sin(x^2)
g = cos(x)


@test_approx_eq_eps f(.1) sin(.1^2) 1000eps()

h = f + g^2
r = roots(h)
rp = roots(differentiate(h))

@test norm(h(r))<1000eps()

@test norm(diff(h)(rp))<100000eps()




## Differentiation and Integration
f = Fun(x->exp(x),[-1.,1.])
fp = differentiate(f)
@test norm(f-fp)<1000eps()

g = cumsum(f)
g = g + f(-1)
@test norm(f-g)<100eps()



x = Fun(identity)
f = exp(x)
g = f/sqrt(1-x^2)


space(f),domain(f)
space(g),domain(g)


f = Fun(x->cospi(5x))
g = abs(f)
space(f)
space(g)


## Solving ODEs

x = Fun(identity,[-1000.,200.])
d = domain(x)
D = Derivative(d)
B = dirichlet(d)
L = D^2 - x
u = [B;L] \ [airyai(d.a);airyai(d.b)]

@test_approx_eq u(0.) airyai(0.)


## Nonlinear BVPs
x=Fun()
u0=0.x

N=u->[u(-1.)-1.,u(1.)+0.5,0.001u''+6*(1-x^2)*u'+u^2-1.]
u=newton(N,u0)

@test norm(N(u)[end]) ≤ 1000eps()

## Periodic Functions

f = Fun(cos,Fourier([-π,π]))
@test norm(differentiate(f) + Fun(sin,Fourier([-π,π])))<100eps()




s = Chebyshev([-π,π])
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
B = periodic(s,0)
uChebyshev = [B;L]\[0.,f]

s = Fourier([-π,π])
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
uFourier = L\f


@ test_approx_eq uChebyshev(0.) uFourier(0.)


## Sampling

f = abs(Fun(sin,[-5,5]))
x = ApproxFun.sample(f,10000)


## PDEs

d = Interval()^2                            # Defines a rectangle

u = [dirichlet(d);lap(d)+100I]\ones(4)      # First four entries of rhs are





d = Interval()^2
u0 = Fun((x,y)->exp(-40(x-.1)^2-40(y+.2)^2),d)
B = dirichlet(d)
D = Derivative(Interval())
L = (0.01D^2-4D)⊗I + I⊗(0.01D^2-3D)
h = 0.002
