using ApproxFun, Base.Test


## Avazzadev et al

# Example 1

x=Fun([0.,1.])
Q=gamma(0.5)*LeftIntegral(0.5)
u=Q\(2/105*sqrt(x)*(105-56x^2+48x^3))
@test norm(u-(x^3-x^2+1))<10eps()


# Example 2

x=Fun([0.,1.])
Q=gamma(0.5)*LeftIntegral(0.5)
u=Q\(exp(x)-1)
@test norm(u-exp(x)*erf(sqrt(x))/sqrt(π)) < 10eps() # 5.0036177384681187e-14


# Example 3
x=Fun([0.,1.])
Q=gamma(1/5)*LeftIntegral(1/5)
u=Q\(x+1)
norm(u.-(1+1.25x)*sin(0.8π)/(π*x^(1/5)))

# Example 4
x=Fun([0.,1.])
Q=gamma(1-1/3)*LeftIntegral(1-1/3)
u=Q\x^(7/6)
norm(u-7*gamma(1/6)/(18*sqrt(π)*gamma(2/3))*sqrt(x))


# Example 5

d=Interval([0.,1.])
x=Fun(d)
f=x+4/3*x^(3/2)
S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
Q=gamma(.5)*LeftIntegral(S,.5)
L=I+Q
u=L\f
norm(u-x)  # 4.967828369398805e-15

# Example 6

d=Interval([0.,1.])
x=Fun(d)
f=x^2+16/15*x^(5/2)
S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
Q=gamma(.5)*LeftIntegral(S,.5)
L=I+Q
u=L\f
norm(u-x^2) # 3.741677490621214e-15

# Example 7

d=Interval([0.,1.])
x=Fun(d)
f=2sqrt(x)
S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
Q=gamma(.5)*LeftIntegral(S,.5)
L=I+Q
u=L\f
plot(u)

norm(1-exp(π*x)*erfc(sqrt(π*x))-u)


# Example 8

d=Interval([0.,1.])
x=Fun(d)
f=1/(x+1)+2*Fun(x->asinh(sqrt(x))/sqrt(1+x),JacobiWeight(.5,0.,d))
S=Legendre(d)⊕JacobiWeight(.5,0.,Jacobi(.5,.5,d))
Q=gamma(.5)*LeftIntegral(S,.5)
L=I+Q
u=L\f
norm((u-1/(x+1)).coefficients)   # 1.2011889731154679e-14
