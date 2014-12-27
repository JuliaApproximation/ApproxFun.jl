using ApproxFun, Base.Test

u0   = TensorFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>LowRankFun)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)[.1,.2]-sin(u0[.1,.2])|>abs < 10e-6


## Rectangle PDE

dx=dy=Interval()
d=dx*dy
x=Fun(identity,dx);y=Fun(identity,dy)

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

G=[real(exp(-1+1.im*y)),
                        real(exp(1+1im*y)),
                        real(exp(x-1im)),
                        real(exp(x+1im)),0.];

A=[dirichlet(d),lap(d)]

u=A\G
@test_approx_eq u[.1,.2] real(exp(.1+.2im))



## Periodic
f=LowRankFun((x,y)->cos(x)*sin(y),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(.1)*sin(.2)

f=LowRankFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=points(f,1),y=points(f,2)]-values(f))<10000eps()

f=TensorFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval()^2)
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=points(f,1),y=points(f,2)]-values(f))<10000eps()

d=PeriodicInterval()^2
f=TensorFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
@test (f.'-f|>coefficients|>norm)< 10eps()



d=PeriodicInterval()^2
f=TensorFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 10000eps()

@test_approx_eq real(f)[.1,.2] f[.1,.2]





## Kron

dx=dy=Interval()
d=dx*dy
x=Fun(identity,dx);y=Fun(identity,dy)

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

G=[real(exp(-1+1.im*y)),
                        real(exp(1+1im*y)),
                        real(exp(x-1im)),
                        real(exp(x+1im)),0.];

A=[dirichlet(d),lap(d)]

S=schurfact(A,40)

uex=A\G

nx=ny=40
K=kron(A,nx,ny)

uex2=K\G

@test (uex-uex2|>coefficients|>norm)<100eps()



# dirichlet bcs

import ApproxFun.ChebyshevDirichlet

S=ChebyshevDirichlet()⊗ChebyshevDirichlet();
A=[dirichlet(S),lap(S)]
nx=ny=20;
KD=kron(A,nx,ny);


#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
x=Fun(identity);y=Fun(identity);
G=[Fun(real(exp(-1+1.im*y)),S[2]),
    Fun(real(exp(1+1im*y)),S[2]),
    Fun(real(exp(x-1im)),S[1]),
                        Fun(real(exp(x+1im)),S[1]),0.];

uD=KD\G;

@test_approx_eq uD[.1,.2] real(exp(.1+.2im))

