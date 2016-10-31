using ApproxFun, Base.Test
    import ApproxFun: testbandedbelowoperator, testbandedoperator, testspace, testtransforms


testspace(Jacobi(2.,.5);haslineintegral=false)

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
testtransforms(JacobiWeight(0.,m,Jacobi(2m+1,0.)))
f=Fun(x->((1-x)/2).^m.*exp(x),JacobiWeight(0.,m,Jacobi(2m+1,0.)))
@test abs(f(.1)-(x->((1-x)/2).^m.*exp(x))(.1))<10eps()


m=10
@time f=Fun(x->besselj(m,m*(1-x)),JacobiWeight(0.,m,Jacobi(2m+1,0.)))
@test_approx_eq f(0.) besselj(m,m)






## Conversion

testtransforms(Jacobi(-0.5,-0.5))
@test norm(Fun(Fun(exp),Jacobi(-.5,-.5))-Fun(exp,Jacobi(-.5,-.5))) < 100eps()

x=Fun(identity)
ri=0.5./(1-x)
@test_approx_eq ((1-x)./2.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) (1-.1)./2*exp(.1)


@test_approx_eq ((1-x)./2.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) (1-.1)./2*exp(.1)


@test_approx_eq (ri.*Fun(exp,JacobiWeight(0.,0.,Jacobi(1.,0.))))(.1) .5/(1-.1)*exp(.1)


## Derivative

D=Derivative(Jacobi(1.,0.,Interval(1.,0.)))
testbandedoperator(D)

S=JacobiWeight(0.,0.,Jacobi(1.,0.,Interval(1.,0.)))
D=Derivative(S)
testbandedoperator(D)

f=Fun(exp,domainspace(D))
@test (D*f-f).coefficients|>norm < eps(100000.)
@test (f'-f).coefficients|>norm < eps(100000.)
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

testtransforms(Legendre([0,2]))
@test_approx_eq sum(Fun(exp,Legendre([0,2]))) sum(Fun(exp,[0,2]))

a=Arc(0.,.1,0.,π/2)
g=Fun(exp,Legendre(a))

@test_approx_eq sum(g) sum(Fun(exp,a))



## Test special derivative

x=Fun()
f=exp(x)*sqrt(1-x^2)
D=Derivative(WeightedJacobi(.5,.5))

testtransforms(WeightedJacobi(.5,.5))
testbandedoperator(D)

@time g=(D*Fun(f,domainspace(D)))
@test_approx_eq f'(0.1) g(0.1)

## Test implementation of conversion between Chebyshev and Jacobi spaces using FastTransforms

f = Fun(x->cospi(1000x))
g = Fun(f,Legendre())
h = Fun(g,Chebyshev())
@test norm(f.coefficients-h.coefficients,Inf) < 100eps()
@time h = Fun(h,Legendre())
@test norm(g.coefficients-h.coefficients,Inf) < 1000eps()




## Legendre conversions
testspace(Ultraspherical(1);haslineintegral=false)
testspace(Ultraspherical(2);haslineintegral=false)
testspace(Ultraspherical(1//2);haslineintegral=false,minpoints=2)  # minpoints is a tempory fix a bug

@test norm(Fun(exp,Ultraspherical(1//2))-Fun(exp,Jacobi(0,0))) < 100eps()

C=Conversion(Jacobi(0,0),Chebyshev())
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Jacobi(0,0))  - Fun(exp)) < 100eps()


C=Conversion(Ultraspherical(1//2),Chebyshev())
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Ultraspherical(1//2))  - Fun(exp)) < 100eps()



C=Conversion(Chebyshev(),Ultraspherical(1//2))
testbandedbelowoperator(C)
@test norm(C*Fun(exp)-Fun(exp,Legendre())) < 100eps()


C=Conversion(Chebyshev(),Jacobi(0,0))
testbandedbelowoperator(C)
@test norm(C*Fun(exp)  - Fun(exp,Jacobi(0,0))) < 100eps()


C=Conversion(Chebyshev(),Jacobi(1,1))
testbandedbelowoperator(C)
@test norm(C*Fun(exp) - Fun(exp,Jacobi(1,1))) < 100eps()


C=Conversion(Ultraspherical(1//2),Ultraspherical(1))
testbandedbelowoperator(C)

λ1 = ApproxFun.order(domainspace(C))
λ2 = ApproxFun.order(rangespace(C))

# test against version that doesn't use lgamma
Cex = Float64[(if j ≥ k && iseven(k-j)
        gamma(λ2)*(k-1+λ2)/(gamma(λ1)*gamma(λ1-λ2))*
            (gamma((j-k)/2+λ1-λ2)/gamma((j-k)/2+1))*
            (gamma((k+j-2)/2+λ1)/gamma((k+j-2)/2+λ2+1))
    else
        0.0
    end) for k=1:20,j=1:20]

@test norm(Cex - C[1:20,1:20]) < 100eps()

@test norm(C*Fun(exp,Ultraspherical(1//2))-Fun(exp,Ultraspherical(1))) < 100eps()

C=Conversion(Jacobi(0,0),Ultraspherical(1))
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Jacobi(0,0))-Fun(exp,Ultraspherical(1))) < 100eps()


C=Conversion(Ultraspherical(1),Jacobi(0,0))
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Ultraspherical(1))-Fun(exp,Jacobi(0,0))) < 100eps()



## Derivative

S=JacobiWeight(1,1,Ultraspherical(1))

f=Fun([1.,2.,3.],S)
@test_approx_eq (Derivative(S,2)*f)(0.1) f''(0.1)


## == tests


@test WeightedJacobi(0.1,0.2) == WeightedJacobi(0.1+eps(),0.2)

# this tests a subspace bug
f=Fun(rand(10),WeightedJacobi(0.1,0.2))  # convert to Legendre expansion

g=(f|(2:ApproxFun.∞))

@test ApproxFun.coefficients(g.coefficients,space(g),ApproxFun.canonicalspace(g))[1] ==0.
@test norm((Fun(g,space(f))|(2:ApproxFun.∞)-g).coefficients) < 10eps()


## Check conversion for non-compatible paramters
S=Jacobi(0.1,1.2)
x=Fun()


p=(S,k)->Fun([zeros(k);1.],S)
n=1;
@test norm(x*p(S,n-1)-(ApproxFun.recα(Float64,S,n)*p(S,n-1) + ApproxFun.recβ(Float64,S,n)*p(S,n))) < 10eps()



# Log with squareroot singularities

a=1.0;b=2.0+im
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
S=space(f)
x=4.0+2im;
@test_approx_eq linesum(f*log(abs(x-z))) 13.740676344264614



# Line sum for legendre

x=Fun(Legendre())
@test_approx_eq sum(x+1) linesum(x+1)

x=Fun(Legendre([2,1]))
@test_approx_eq sum(x+1) -linesum(x+1)

x=Fun(Interval(1,1+im))
@test_approx_eq sum(x+1) im*linesum(x+1)

x=Fun(Legendre(Interval(1,1+im)))
@test_approx_eq sum(x+1) im*linesum(x+1)

x=Fun(Legendre(Interval(im,1)))
sum(x+1), linesum(x+1)
