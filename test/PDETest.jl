using ApproxFun, Compat, Base.Test

## Rectangle PDE

dx=dy=Interval()
d=dx*dy
x=Fun(identity,dx);y=Fun(identity,dy)

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

G=[real(exp(-1+1.0im*y));
                        real(exp(1+1im*y));
                        real(exp(x-1im));
                        real(exp(x+1im));0.];

A=[dirichlet(d);lap(d)]
u=A\G
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))


A=[dirichlet(d);lap(d)+0.0I]
u=A\G
@test_approx_eq_eps u(.1,.2) real(exp(0.1+0.2im)) 1E-11


println("    Poisson tests")

## Poisson

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2))  #default is [-1,1]^2
d=domain(f)
OS=S=schurfact([dirichlet(d);lap(d)],10)
u=OS\[zeros(∂(d));f]
@test_approx_eq u(.1,.2) -0.042393137972085826



d=PeriodicInterval()^2
f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 1000000eps()


println("    Kron tests")

@static if is_apple()
    ## Kron

    dx=dy=Interval()
    d=dx*dy
    x=Fun(identity,dx);y=Fun(identity,dy)

    #dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]

    G=[real(exp(-1+1.0im*y));
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
    G=[Fun(real(exp(-1+1.0im*y)),S[2]);
        Fun(real(exp(1+1im*y)),S[2]);
        Fun(real(exp(x-1im)),S[1]);
                            Fun(real(exp(x+1im)),S[1]);0.];

    uD=KD\G;

    @test_approx_eq uD(.1,.2) real(exp(.1+.2im))



    # fourth order
    dx=dy=Interval()
    d=dx*dy
    Dx=Derivative(dx);Dy=Derivative(dy)
    L=Dx^4⊗I+2*Dx^2⊗Dy^2+I⊗Dy^4

    K=kronfact([dirichlet(d);
         neumann(d);
         L],100,100)

    x=Fun(identity,dx);y=Fun(identity,dy)

    G=[real(exp(-1+1.0im*y));
                    real(exp(1+1im*y));
                    real(exp(x-1im));
                    real(exp(x+1im));
                    real(exp(-1+1.0im*y));
                    real(exp(1+1im*y));
                    -imag(exp(x-1im));
                    -imag(exp(x+1im))
       ]
    u=K\G
    @test_approx_eq u(.1,.2) real(exp(.1+.2im))


    # mixed

    K=kronfact([(ldirichlet(dx)+lneumann(dx))⊗I;
            (rdirichlet(dx)+rneumann(dx))⊗I;
            I⊗(ldirichlet(dy)+lneumann(dy));
            I⊗(rdirichlet(dy)+rneumann(dy));
            (ldirichlet(dx)-lneumann(dx))⊗I;
            (rdirichlet(dx)-rneumann(dx))⊗I;
            I⊗(ldirichlet(dy)-lneumann(dy));
            I⊗(rdirichlet(dy)-rneumann(dy));
             L],100,100)
    G=[2real(exp(-1+1.0im*y));
                    2real(exp(1+1im*y));
                    real(exp(x-1im))-imag(exp(x-1im));
                    real(exp(x+1im))-imag(exp(x+1im));
                    0;
                    0;
                    real(exp(x-1im))+imag(exp(x-1im));
                    real(exp(x+1im))+imag(exp(x+1im))
       ]
    u=K\G

    @test_approx_eq u(.1,.2) real(exp(.1+.2im))
end



## Test periodic x interval

d=PeriodicInterval()*Interval()
g=Fun(z->real(cos(z)),∂(d))  # boundary data
u=[dirichlet(d);lap(d)]\g

@test_approx_eq u(.1,.2) real(cos(.1+.2im))



dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
u=[I⊗ldirichlet(dt);Dt+Dθ]\Fun(θ->exp(-20θ^2),dθ)


@test_approx_eq u(.1,.2) exp(-20(0.1-0.2)^2)


println("   Domain Decompositon tests")

## Domain Decomposition
d=Interval(0,1)^2
n,m=20,80
A=discretize([dirichlet(d);lap(d)],n)
∂d=∂(d)
g=Fun(z->real(exp(z)),∂d)
f=[Fun([zeros(k-1);1.0],∂d) for k=1:m].'
@time U=A\f
@test_approx_eq dot(real(g.coefficients),U[1:ncoefficients(g)])(.1,.2) real(exp(.1+.2im))



Rectangle(a,b,c,d)=Interval(a,b)*Interval(c,d)
Γ=Rectangle(0,1,0,1)∪Rectangle(1,2,0,1)
Fun(identity,∂(Γ))|>values



## Small diffusoion
println("   Small diffusion tests")

dx=Interval();dt=Interval(0,1.)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x=Fun(identity,dx)
B=0.0
C=0.0
V=B+C*x
ε=0.001
f=Fun(x->exp(-20x^2),dx)
@time u=[timedirichlet(d);Dt-ε*Dx^2-V*Dx]\f
@test_approx_eq u(.1,.2) 0.8148207991358946
B=0.1
C=0.2
V=B+C*x
@time u=[timedirichlet(d);Dt-ε*Dx^2-V*Dx]\f
@test_approx_eq u(.1,.2) 0.7311625132209619


## Schrodinger

dx=Interval(0.,1.);dt=Interval(0.0,.1)
d=dx*dt

V=Fun(x->x^2,dx)

Dt=Derivative(d,[0,1]);Dx=Derivative(d,[1,0])

ϵ=1.
u0=Fun(x->exp(-100*(x-.5)^2)*exp(-1./(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
L=ϵ*Dt+(.5im*ϵ^2*Dx^2)
ny=200;u=pdesolve([timedirichlet(d);L],u0,ny)
@test_approx_eq_eps u(.2,.1) (0.2937741918470843 + 0.22130344715160255im )  0.000001



## Periodic

d=PeriodicInterval()^2
f=Fun((θ,ϕ)->exp(-10(sin(θ/2)^2+sin(ϕ/2)^2)),d)
A=lap(d)+.1I
u=A\f
@test_approx_eq u(.1,.2) u(.2,.1)


d=PeriodicInterval()*Interval()
g=Fun(z->real(cos(z)),∂(d))  # boundary data
u=[dirichlet(d);lap(d)]\g
@test_approx_eq u(.1,.2) real(cos(.1+.2im))



dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
u=[I⊗ldirichlet(dt);Dt+Dθ]\Fun(θ->exp(-20θ^2),dθ)




# Check bug in cache
CO=cache(ldirichlet(dt))
ApproxFun.resizedata!(CO,:,2)
ApproxFun.resizedata!(CO,:,4)
@test_approx_eq CO*Fun(exp,dt) 1.0

dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dt*dθ
Dt=Derivative(d,[1,0]);Dθ=Derivative(d,[0,1])
A=[ldirichlet(dt)⊗I;Dt+Dθ]
f=Fun(θ->exp(-20θ^2),dθ)
ut=A\f

@test_approx_eq u(.1,.2) ut(.2,.1)




# Beam

dθ=PeriodicInterval(0.0,1.0);dt=Interval(0,0.03)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1]);

B=[I⊗ldirichlet(dt),I⊗lneumann(dt)]
u=pdesolve([B;Dt^2+Dθ^4],Fun(θ->exp(-200(θ-.5).^2),dθ),200)

@test_approx_eq_eps u(.1,.01) -0.2479768394633227  1E-8 #empirical



## Rectangle PDEs

# Screened Poisson

d=Interval()^2
u=[neumann(d);lap(d)-100.0I]\ones(∂(d))
@test_approx_eq u(.1,.9) 0.03679861429138079

# PiecewisePDE

a=Fun([1,0.5,1],[-1.,0.,0.5,1.])
s=space(a)
dt=Interval(0,2.)
Dx=Derivative(s);Dt=Derivative(dt)
Bx=[ldirichlet(s);continuity(s,0)]


# test resize bug
CO=cache(Bx[2])
ApproxFun.resizedata!(CO,:,2)
ApproxFun.resizedata!(CO,:,4)
@test_approx_eq CO.data*collect(1:4) [3.,-1.]


u=pdesolve([I⊗ldirichlet(dt);Bx⊗I;I⊗Dt+(a*Dx)⊗I],Any[Fun(x->exp(-20(x+0.5)^2),s)],200)
@test_approx_eq_eps u(-.1,.2) exp(-20(-.2-.1+0.5)^2) 0.00001



## Test error


dx=Interval();dt=Interval(0,2.)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x=Fun(identity,dx)
u=[I⊗ldirichlet(dt);Dt+x*Dx]\Fun(x->exp(-20x^2),dx)

@test_approx_eq u(0.1,0.2) 0.8745340845783758  # empirical


dθ=PeriodicInterval();dt=Interval(0,10.)
d=dθ*dt
ε=.01
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])

# Parentheses are a hack to get rank 2 PDE
u=[I⊗ldirichlet(dt);Dt-ε*Dθ^2-Dθ]\Fun(θ->exp(-20θ^2),dθ)

@test_approx_eq_eps u(0.1,0.2) 0.1967278179230314 1000eps()
