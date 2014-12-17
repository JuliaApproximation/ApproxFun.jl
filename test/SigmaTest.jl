using ApproxFun, Base.Test

#The first test checks the solution of the integral equation
# u(x) + \int_{-1}^{+1} \frac{e^{y} u(y)}{\sqrt{1-y^2}} dy = f
# on the interval [-1,1].

x=Fun(identity)
w=1/sqrt(1-x^2)
d=domain(x)

S=Σ(d)

@test domainspace(S) == JacobiWeightSpace{ChebyshevSpace}(-0.5,-0.5,ChebyshevSpace())
@test rangespace(S) == ChebyshevSpace()

L=I+S[exp(x)*w]
usol=sin(2x)
f=L*usol
u=L\f
@test norm(u-usol) <= 10eps()


#The second test checks the solution of the integro-differential equation
# u'(x) + x u(x) + \int_{-2}^{+2} sin(y-x) u(y) \sqrt{4-y^2} dy = f
# on the interval [-2,2], with u(-2) = 1.

x=Fun(identity,[-2.,2.])
w=sqrt(4-x^2)
d=domain(x)

D=Derivative(d)
B=ldirichlet(d)
S=Σ(.5,.5,d)

@test domainspace(S) == JacobiWeightSpace{UltrasphericalSpace{1}}(.5,.5,UltrasphericalSpace{1}(d))
@test rangespace(S) == UltrasphericalSpace{1}(d)

K=LowRankFun((x,y)->sin(y-x)*w[y],rangespace(S),domainspace(S))

L=D+x+S[K]
usol=cospi(20x)
f=L*usol
u=[B,L]\[1.,f]
@test norm(u-usol) <= 10eps()