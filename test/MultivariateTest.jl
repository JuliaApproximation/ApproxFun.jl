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
    K=kron(A,nx,ny)

    uex2=K\G

    @test (uex-uex2|>coefficients|>norm)<10000eps()



    # dirichlet bcs

    import ApproxFun.ChebyshevDirichlet

    S=ChebyshevDirichlet()⊗ChebyshevDirichlet();
    A=[dirichlet(S);lap(S)]
    nx=ny=20;
    KD=kron(A,nx,ny);


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




