using ApproxFun, Base.Test



## Try constructor variants

ff=(x,y)->exp(-10(x+.2)^2-20(y-.1)^2)*cos(x*y)
gg=x->exp(-10(x[1]+.2)^2-20(x[1]-.1)^2)*cos(x[1]*x[2])
f=Fun(ff,Interval()^2,10000)
@test_approx_eq f[0.,0.] ff(0.,0.)

f=Fun(gg,Interval()^2,10000)
@test_approx_eq f[0.,0.] ff(0.,0.)

f=Fun(ff,Interval()^2)
@test_approx_eq f[0.,0.] ff(0.,0.)
f=Fun(gg,Interval()^2)
@test_approx_eq f[0.,0.] ff(0.,0.)


f=Fun(ff)
@test_approx_eq f[0.,0.] ff(0.,0.)
f=Fun(gg)
@test_approx_eq f[0.,0.] ff(0.,0.)



## ProductFun
u0   = ProductFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>LowRankFun)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)[.1,.2]-sin(u0[.1,.2])|>abs < 10e-4

## LowRankFun

F = LowRankFun((x,y)->besselj0(10(y-x)),Chebyshev(),Chebyshev())

@test_approx_eq F[.123,.456] besselj0(10(.456-.123))

G = LowRankFun((x,y)->besselj0(10(y-x));method=:Cholesky)

@test_approx_eq G[.357,.246] besselj0(10(.246-.357))

F = LowRankFun((x,y)->hankelh1(0,10abs(y-x)),Chebyshev([1.0,2.0]),Chebyshev([1.0im,2.0im]))

@test_approx_eq F[1.5,1.5im] hankelh1(0,10abs(1.5im-1.5))

## Rectangle PDE

dx=dy=Interval()
d=dx*dy
x=Fun(identity,dx);y=Fun(identity,dy)

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

G=[real(exp(-1+1.im*y));
                        real(exp(1+1im*y));
                        real(exp(x-1im));
                        real(exp(x+1im));0.];

A=[dirichlet(d);lap(d)]
u=A\G
@test_approx_eq u[.1,.2] real(exp(.1+.2im))


## Poisson

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2))  #default is [-1,1]^2
d=domain(f)
OS=S=schurfact([dirichlet(d);lap(d)],10)
u=OS\[zeros(∂(d));f]
@test_approx_eq u[.1,.2] -0.042393137972085826



## Periodic
f=LowRankFun((x,y)->cos(x)*sin(y),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(.1)*sin(.2)

f=LowRankFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=ApproxFun.vecpoints(f,1),y=ApproxFun.vecpoints(f,2)]-values(f))<10000eps()

f=ProductFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval()^2)
@test_approx_eq f[.1,.2] cos(cos(.1)+sin(.2))
x,y=points(f)
@test norm(Float64[cos(cos(x[k,j])+sin(y[k,j])) for k=1:size(f,1),j=1:size(f,2)]-values(f))<10000eps()

d=PeriodicInterval()^2
f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
@test (f.'-f|>coefficients|>norm)< 1000eps()



d=PeriodicInterval()^2
f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 1000000eps()

@test_approx_eq real(f)[.1,.2] f[.1,.2]



if OS_NAME==:Darwin
    ## Kron

    dx=dy=Interval()
    d=dx*dy
    x=Fun(identity,dx);y=Fun(identity,dy)

    #dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

    G=[real(exp(-1+1.im*y));
                            real(exp(1+1im*y));
                            real(exp(x-1im));
                            real(exp(x+1im));0.];

    A=[dirichlet(d);lap(d)]

    S=schurfact(A,40)

    uex=A\G

    nx=ny=40
    K=kronfact(A,nx,ny)

    uex2=K\G

    @test (uex-uex2|>coefficients|>norm)<10000eps()



    # dirichlet bcs

    import ApproxFun.ChebyshevDirichlet

    S=ChebyshevDirichlet()⊗ChebyshevDirichlet();
    A=[dirichlet(S);lap(S)]
    nx=ny=20;
    KD=kronfact(A,nx,ny);


    #dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
    x=Fun(identity);y=Fun(identity);
    G=[Fun(real(exp(-1+1.im*y)),S[2]);
        Fun(real(exp(1+1im*y)),S[2]);
        Fun(real(exp(x-1im)),S[1]);
                            Fun(real(exp(x+1im)),S[1]);0.];

    uD=KD\G;

    @test_approx_eq uD[.1,.2] real(exp(.1+.2im))
end



## Test periodic x interval

d=PeriodicInterval()*Interval()
g=Fun(z->real(cos(z)),∂(d))  # boundary data
u=[dirichlet(d);lap(d)]\g

@test_approx_eq u[.1,.2] real(cos(.1+.2im))



dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dθ*dt
Dθ=diff(d,1);Dt=diff(d,2)
u=[I⊗ldirichlet(dt);Dt+Dθ]\Fun(θ->exp(-20θ^2),dθ)


@test_approx_eq u[.1,.2] exp(-20(0.1-0.2)^2)




## Functional*Fun

d=Interval()
B=ldirichlet(d)
f=ProductFun((x,y)->cos(cos(x)*sin(y)),d^2)
@test norm(B*f-Fun(y->cos(cos(-1)*sin(y)),d))<20000eps()
@test norm(f*B-Fun(x->cos(cos(x)*sin(-1)),d))<20000eps()



## Domain Decomposition



d=Interval(0,1)^2
A=discretize([dirichlet(d);lap(d)],20)
∂d=∂(d)
g=Fun(z->real(exp(z)),∂d)
f=[Fun([zeros(k-1);1.0],∂d) for k=1:80].'
U=A\f
@test_approx_eq dot(real(g.coefficients),U[1:length(g)])[.1,.2] real(exp(.1+.2im))




## Disk

# This disables this test since it depends on some lib
if OS_NAME == :Darwin
    # Laplace
    d=Disk()
    u=[dirichlet(d),lap(d)]\Fun(z->real(exp(z)),Circle())
    @test_approx_eq u[.1,.2] real(exp(.1+.2im))

    # remaining numbers determined numerically, may be
    # inaccurate

    # Poisson
    f=Fun((x,y)->exp(-10(x+.2).^2-20(y-.1).^2),d)
    u=[dirichlet(d),lap(d)]\[0.,f]
    @test_approx_eq u[.1,.2] -0.039860694987858845

    #Helmholtz
    u=[dirichlet(d),lap(d)+100I]\1.0
    @test_approx_eq u[.1,.2] -0.3675973169667076
    u=[neumann(d),lap(d)+100I]\1.0
    @test_approx_eq u[.1,.2] -0.20795862954551195

    # Screened Poisson
    u=[neumann(d),lap(d)-100.0I]\1.0
    @test_approx_eq u[.1,.9] 0.04313812031635443

    # Lap^2
    u=[dirichlet(d),neumann(d),lap(d)^2]\Fun(z->real(exp(z)),Circle())
    @test_approx_eq u[.1,.2] 1.1137317420521624
end




## Small diffusoion

dx=Interval();dt=Interval(0,1.)
d=dx*dt
Dx=diff(d,1);Dt=diff(d,2)
x=Fun(identity,dx)
B=0.0
C=0.0
V=B+C*x
ε=0.001
f=Fun(x->exp(-20x^2),dx)
u=[timedirichlet(d);Dt-ε*Dx^2-V*Dx]\f


@test_approx_eq u[.1,.2] 0.8148207991358946
B=0.1
C=0.2
V=B+C*x
u=[timedirichlet(d);Dt-ε*Dx^2-V*Dx]\f
@test_approx_eq u[.1,.2] 0.7311625132209619


## Schrodinger

dx=Interval(0.,1.);dt=Interval(0.0,.1)
d=dx*dt

V=Fun(x->x^2,dx)

Dt=diff(d,2);Dx=diff(d,1)

ϵ=1.
u0=Fun(x->exp(-100*(x-.5)^2)*exp(-1./(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
L=ϵ*Dt+(.5im*ϵ^2*Dx^2)
ny=200;u=pdesolve([timedirichlet(d);L],u0,ny)
@test_approx_eq_eps u[.2,.1] (0.2937741918470843 + 0.22130344715160255im )  0.000001



## Periodic

d=PeriodicInterval()^2
f=Fun((θ,ϕ)->exp(-10(sin(θ/2)^2+sin(ϕ/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test_approx_eq u[.1,.2] u[.2,.1]


d=PeriodicInterval()*Interval()
g=Fun(z->real(cos(z)),∂(d))  # boundary data
u=[dirichlet(d);lap(d)]\g
@test_approx_eq u[.1,.2] real(cos(.1+.2im))



dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dθ*dt
Dθ=Derivative(d,1);Dt=Derivative(d,2)
u=[I⊗ldirichlet(dt);Dt+Dθ]\Fun(θ->exp(-20θ^2),dθ)

A=[ldirichlet(dt)⊗I;(Dt+Dθ).']
f=Fun(θ->exp(-20θ^2),dθ)
ut=A\f

@test_approx_eq u[.1,.2] ut[.2,.1]




# Beam

dθ=PeriodicInterval(0.0,1.0);dt=Interval(0,0.03)
d=dθ*dt
Dθ=Derivative(d,1);Dt=Derivative(d,2);

B=[I⊗ldirichlet(dt),I⊗lneumann(dt)]
u=pdesolve([B;Dt^2+Dθ^4],Fun(θ->exp(-200(θ-.5).^2),dθ),200)

@test_approx_eq_eps u[.1,.01] -0.2479768394633227  1E-8 #empirical



## Rectangle PDEs

# Screened Poisson

d=Interval()^2
u=[neumann(d);lap(d)-100.0I]\ones(∂(d))
@test_approx_eq u[.1,.9] 0.03679861429138079

# PiecewisePDE

a=Fun([1,0.5,1],[-1.,0.,0.5,1.])
s=space(a)
dt=Interval(0,2.)
Dx=Derivative(s);Dt=Derivative(dt)
Bx=[ldirichlet(s);continuity(s,0)]
u=pdesolve([I⊗ldirichlet(dt);Bx⊗I;I⊗Dt+(a*Dx)⊗I],Any[Fun(x->exp(-20(x+0.5)^2),s)],200)
       #.'
@test_approx_eq_eps u[-.1,.2] exp(-20(-.2-.1+0.5)^2) 0.00001

