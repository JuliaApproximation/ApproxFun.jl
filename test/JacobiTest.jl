using ApproxFun, Base.Test


f=Fun(exp,Jacobi(2.,.5))
@test_approx_eq f(.1) exp(.1)

f=Fun(x->cos(100x),Jacobi(2.124,.5),500)
@test_approx_eq f(.1) cos(100*.1)


sp=Jacobi(2.124,.5)
f=Fun(exp,sp)
sp2=Jacobi(2.124,1.5)
f2=Fun(exp,sp2)
sp3=Jacobi(3.124,1.5)
f3=Fun(exp,sp3)
sp4=Jacobi(4.124,2.5)
f4=Fun(exp,sp4)
@test norm((Fun(f,sp2)-f2).coefficients)<10eps()
@test norm((Fun(f,sp3)-f3).coefficients)<10eps()
@test norm((Fun(f,sp4)-f4).coefficients)<20eps()


m=20
f=Fun(x->((1-x)/2).^m.*exp(x),JacobiWeight(0.,m,Jacobi(2m+1,0.)))
@test abs(f(.1)-(x->((1-x)/2).^m.*exp(x))(.1))<10eps()


m=10
f=Fun(x->besselj(m,m*(1-x)),JacobiWeight(0.,m,Jacobi(2m+1,0.)))
@test_approx_eq f(0.) besselj(m,m)






## Conversion

@test norm(Fun(Fun(exp),Jacobi(-.5,-.5))-Fun(exp,Jacobi(-.5,-.5))) < 100eps()

x=Fun(identity)
ri=0.5./(1-x)
@test_approx_eq ((1-x)./2.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) (1-.1)./2*exp(.1)


@test_approx_eq ((1-x)./2.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) (1-.1)./2*exp(.1)


@test_approx_eq (ri.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) .5/(1-.1)*exp(.1)


## Derivative

S=JacobiWeight(0.,0.,Jacobi(1.,0.,Interval(1.,0.)))
D=Derivative(S)
f=Fun(exp,domainspace(D))
@test (D*f-f).coefficients|>norm < eps(100000.)
@test (diff(f)-f).coefficients|>norm < eps(100000.)
@test (D^2*f-f).coefficients|>norm < eps(100000000.)
@test (D*(D*f)-f).coefficients|>norm < eps(100000000.)





### Jacobi multiplication

x=Fun(identity,Jacobi(0.,0.))
f=Fun(exp,Jacobi(0.,0.))

@test_approx_eq (x*f)(.1) .1exp(.1)

x=Fun(identity,Jacobi(0.123,12.324))
f=Fun(exp,Jacobi(0.,0.))

@test_approx_eq (x*f)(.1) .1exp(.1)


x=Fun(identity,Jacobi(0.123,12.324))
f=Fun(exp,Jacobi(0.213,0.590))

@test_approx_eq (x*f)(.1) .1exp(.1)

g=Fun(cos,Jacobi(0.123,12.324))
f=Fun(exp,Jacobi(0.213,0.590))

@test_approx_eq (g*f)(.1) cos(.1)*exp(.1)


## Jacobi integrate and sum

@test_approx_eq sum(Fun(exp,Legendre([0,2]))) sum(Fun(exp,[0,2]))

a=Arc(0.,.1,0.,π/2)
g=Fun(exp,MappedSpace(a,Legendre()))

@test_approx_eq sum(g) sum(Fun(exp,a))
