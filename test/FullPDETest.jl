using ApproxFun, Compat, Base.Test
    import Compat: view
    import ApproxFun: resizedata!, CachedOperator, RaggedMatrix, testbandedblockbandedoperator,
                        testbandedblockoperator, linsolve_coefficients
## Check operators

## Rectangle PDEs

println("    Rectangle tests")


S=JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ=Laplacian(S)

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),rangespace(Δ))  #default is [-1,1]^2
@time v=linsolve(Δ,f;tolerance=1E-14)
@test norm((Δ*v-f).coefficients)<1E-14



# Screened Poisson

dx=dy=Interval()
d=dx*dy
g=Fun((x,y)->exp(x)*cos(y),∂(d))


testbandedblockoperator(Dirichlet(d))


testbandedblockbandedoperator(Laplacian(d))
A=[Dirichlet(d);Laplacian(d)]

@time u=A\[g,0.]
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))



d=Interval()^2
@time u=linsolve([neumann(d);Laplacian(d)-100.0I],[ones(4);0.];tolerance=1E-12)
@test_approx_eq u(.1,.9) 0.03679861429138079



## Test error


dx=Interval();dt=Interval(0,2.)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x,y=Fun(identity,d)
@time u=linsolve([I⊗ldirichlet(dt);Dt+x*Dx],[Fun(x->exp(-20x^2),dx);0.];tolerance=1E-12)

@test_approx_eq u(0.1,0.2) 0.8745340845783758  # empirical




println("    Bilaplacian Tests")

dx=dy=Interval()
d=dx*dy
Dx=Derivative(dx);Dy=Derivative(dy)
L=Dx^4⊗I+2*Dx^2⊗Dy^2+I⊗Dy^4

testbandedblockbandedoperator(L)

A=[dirichlet(dx)⊗eye(dy);
        eye(dx)⊗dirichlet(dy);
        neumann(dx)⊗eye(dy);
        eye(dx)⊗neumann(dy);
         L]


# Checks bug in constructor
f=Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]),22)
@test_approx_eq f(-1.,0.1) real(exp(-1.+0.1im))
f=Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]))
@test_approx_eq f(-1.,0.1) real(exp(-1.+0.1im))

F=[Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[2]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[3]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[4]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[5]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[6]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[7]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[8]));
    0]

@time u=linsolve(A,F;tolerance=1E-10)

@test_approx_eq u(0.1,0.2)  exp(0.1)*cos(0.2)

dx=dy=Interval()
d=dx*dy
Dx=Derivative(dx);Dy=Derivative(dy)
L=Dx^4⊗I+2*Dx^2⊗Dy^2+I⊗Dy^4

testbandedblockbandedoperator(L)

A=[(ldirichlet(dx)+lneumann(dx))⊗eye(dy);
        (rdirichlet(dx)+rneumann(dx))⊗eye(dy);
        eye(dx)⊗(ldirichlet(dy)+lneumann(dy));
        eye(dx)⊗(rdirichlet(dy)+rneumann(dy));
        (ldirichlet(dx)-lneumann(dx))⊗eye(dy);
        (rdirichlet(dx)-rneumann(dx))⊗eye(dy);
        eye(dx)⊗(ldirichlet(dy)-lneumann(dy));
        eye(dx)⊗(rdirichlet(dy)-rneumann(dy));
         L]


u=linsolve(A,[ones(8);0];tolerance=1E-5)
@test_approx_eq u(0.1,0.2) 1.0



F=[2Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]));
    2Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[2]));
    Fun((x,y)->real(exp(x+1.0im*y))-imag(exp(x+1.0im*y)),rangespace(A[3]));
    Fun((x,y)->real(exp(x+1.0im*y))-imag(exp(x+1.0im*y)),rangespace(A[4]));
    0;
    0;
    Fun((x,y)->real(exp(x+1.0im*y))+imag(exp(x+1.0im*y)),rangespace(A[7]));
    Fun((x,y)->real(exp(x+1.0im*y))+imag(exp(x+1.0im*y)),rangespace(A[8]));
    0]

u=linsolve(A,F;tolerance=1E-10)

@test_approx_eq u(0.1,0.2)  exp(0.1)*cos(0.2)


println("    Operator resize tests")

S=ChebyshevDirichlet()^2
B=Dirichlet(S)
f = Fun((x,y)->exp(x)*sin(y),S)
@test norm((Fun((x,y)->exp(x)*sin(y),∂(domain(S))) - B*f).coefficients) < 100eps()


S=JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ=Laplacian(S)

@test_approx_eq cache(Δ)[1:100,1:100]  Δ[1:100,1:100]
@test_approx_eq cache(Δ;padding=true)[1:100,1:100]  Δ[1:100,1:100]

@test_approx_eq cache(Δ)[5:100,7:100]  Δ[5:100,7:100]
@test_approx_eq cache(Δ;padding=true)[5:100,7:100]  Δ[5:100,7:100]

# Check that QR is growing correctly
for col in (1,2,3,10,11,40)
    QR=qrfact(Δ)
    resizedata!(QR.R,:,col+100)
    resizedata!(QR,:,col)
    QR2=qrfact!(CachedOperator(RaggedMatrix,Δ;padding=true))
    resizedata!(QR2.R,:,col+100)
    resizedata!(QR2,:,col)
    n=min(size(QR.H,1),size(QR2.H,1))
    @test_approx_eq QR.H[1:n,1:col] QR2.H[1:n,1:col]
    @test_approx_eq QR.R[1:col,1:col] QR2.R[1:col,1:col]
    @test_approx_eq QR.R[1:col+10,1:col+10] QR2.R[1:col+10,1:col+10]
end

QR=qrfact(Δ)
QR2=qrfact!(CachedOperator(RaggedMatrix,Δ;padding=true))
for col in (80,200)
    resizedata!(QR,:,col)
    resizedata!(QR2,:,col)
    n=min(size(QR.H,1),size(QR2.H,1))
    @test_approx_eq QR.H[1:n,1:col] QR2.H[1:n,1:col]
    @test_approx_eq QR.R[1:col,1:col] QR2.R[1:col,1:col]
    @test_approx_eq QR.R[1:col+10,1:col+10] QR2.R[1:col+10,1:col+10]
end

# this checks a bug
QR=qrfact(Δ)
resizedata!(QR,:,548)
resizedata!(QR,:,430)


u=Fun((x,y)->sin(π*x)*sin(π*y),S)
f=-2π^2*u


QR=qrfact(Δ)
v=QR\f
@test norm((u-v).coefficients)<100eps()

v=Δ\f
@test norm((u-v).coefficients)<100eps()


f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),rangespace(Δ))  #default is [-1,1]^2
@time v=linsolve(Δ,f;tolerance=1E-14)
@test norm((Δ*v-f).coefficients)<1E-14

KO=Δ.op.ops[1].ops[1].op

M=ApproxFun.BandedBlockBandedMatrix(view(KO,1:4,1:4))
@test norm(ApproxFun.BandedBlockBandedMatrix(view(KO,1:4,2:4))-M[:,2:4]) < 10eps()
@test norm(ApproxFun.BandedBlockBandedMatrix(view(KO,1:4,3:4))-M[:,3:4]) < 10eps()

M=ApproxFun.BandedBlockBandedMatrix(view(KO,1:112,1:112))
@test norm(ApproxFun.BandedBlockBandedMatrix(view(KO,1:112,112:112))-M[:,112]) < 10eps()


M=ApproxFun.BandedBlockBandedMatrix(view(Δ,1:4,1:4))
@test norm(ApproxFun.BandedBlockBandedMatrix(view(Δ,1:4,2:4))-M[:,2:4]) < 10eps()
@test norm(ApproxFun.BandedBlockBandedMatrix(view(Δ,1:4,3:4))-M[:,3:4]) < 10eps()

M=ApproxFun.BandedBlockBandedMatrix(view(Δ,1:112,1:112))
@test norm(ApproxFun.BandedBlockBandedMatrix(view(Δ,1:112,112:112))-M[:,112]) < 10eps()



## Rectangle PDE
dx=dy=Interval()
d=dx*dy
g=Fun((x,y)->exp(x)*cos(y),∂(d))

A=[Dirichlet(d);Laplacian(d)]
let Ai=ApproxFun.interlace(A),co=cache(RaggedMatrix,Ai)
    ApproxFun.resizedata!(co,:,100)
    ApproxFun.resizedata!(co,:,200)
    @test norm(Ai[1:200,1:200]-co[1:200,1:200]) == 0
end


u=A\[g,0.]
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))

A=[Dirichlet(d);Laplacian(d)+0.0I]
u=A\[g,0.]
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))



# Check resizing

d=Interval()^2
A=ApproxFun.interlace([Dirichlet(d);Laplacian()+100I])
QR = qrfact(A)
@time ApproxFun.resizedata!(QR.R,:,2000)
@test norm(QR.R.data[1:200,1:200] - A[1:200,1:200]) ==0

@time ApproxFun.resizedata!(QR,:,200)
j=56
v=QR.R.op[1:100,j]
@test norm(linsolve_coefficients(QR[:Q],v;maxlength=300)[j+1:end]) < 100eps()

j=195
v=QR.R.op[1:ApproxFun.colstop(QR.R.op,j),j]
@test norm(linsolve_coefficients(QR[:Q],v;maxlength=1000)[j+1:end]) < 100eps()


j=300
v=QR.R.op[1:ApproxFun.colstop(QR.R.op,j),j]
@test norm(linsolve_coefficients(QR[:Q],v;maxlength=1000)[j+1:end]) < j*20eps()

@test ApproxFun.colstop(QR.R.op,195)-194 == ApproxFun.colstop(QR.H,195)


QR1 = qrfact(A)
@time ApproxFun.resizedata!(QR1.R,:,1000)
QR2 = qrfact([Dirichlet(d);Laplacian()+100I])
@time ApproxFun.resizedata!(QR2.R,:,500)
n=450;QR1.R.data[1:n,1:n]-QR2.R.data[1:n,1:n]|>norm
@time ApproxFun.resizedata!(QR2.R,:,1000)
N=450;QR1.R.data[1:N,1:N]-QR2.R.data[1:N,1:N]|>norm
N=1000;QR1.R.data[1:N,1:N]-QR2.R.data[1:N,1:N]|>norm

QR1 = qrfact(A)
@time ApproxFun.resizedata!(QR1,:,1000)
QR2 = qrfact([Dirichlet(d);Laplacian()+100I])
@time ApproxFun.resizedata!(QR2,:,500)
@time ApproxFun.resizedata!(QR2,:,1000)

@test norm(QR1.H[1:225,1:1000]-QR2.H[1:225,1:1000]) ≤ 10eps()

QR1 = qrfact(A)
@time ApproxFun.resizedata!(QR1,:,5000)
@time u=linsolve(QR1,[ones(∂(d));0.];tolerance=1E-7)

@test norm((Dirichlet(d)*u-ones(∂(d))).coefficients) < 1E-7
@test norm((A*u-Fun([ones(∂(d));0.])).coefficients) < 1E-7
@test norm(((A*u)[2]-(Laplacian()+100I)*u).coefficients) < 1E-10
@test norm((Laplacian()*u+100*u - (A*u)[2]).coefficients) < 1E-10
@time v=linsolve(A,[ones(∂(d));0.];tolerance=1E-7)
@test norm((u-v).coefficients) < 100eps()

@test_approx_eq u(0.1,1.) 1.0
@test_approx_eq u(0.1,-1.) 1.0
@test_approx_eq u(1.,0.1) 1.0
@test_approx_eq u(-1.,0.1) 1.0

S=ChebyshevDirichlet()^2
ff=(x,y)->exp(x)*cos(y)
u=Fun(ff,S)

for KO in [eye(S[1])⊗rdirichlet(S[1]),rdirichlet(S[1])⊗eye(S[2])]
    @test norm((KO*u-Fun(ff,rangespace(KO))).coefficients) ≤ 1E-10
end

B=[dirichlet(S[1])⊗eye(S[2]);
   eye(S[1])⊗dirichlet(S[2]);
   Laplacian()]


u=linsolve(B,[ones(4);0];tolerance=1E-14)
@test norm((u-Fun(S,[1.])).coefficients)<10eps()

g=map(sp->Fun(ff,sp),map(rangespace,B[1:4]))

u=linsolve(B,[g;0];tolerance=1E-10)
@test_approx_eq u(0.1,0.2) ff(0.1,0.2)



println("    Poisson tests")


## Poisson

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),Interval()^2,500)  #default is [-1,1]^2
d=domain(f)
A=[Dirichlet(d);Laplacian(d)]
@time  u=linsolve(A,[zeros(∂(d));f];tolerance=1E-7)
@test_approx_eq_eps u(.1,.2) -0.04251891975068446 1E-5


println("    Periodic Poisson tests")



d=PeriodicInterval()^2
f=Fun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=Laplacian(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 1000000eps()





# fourth order

println("    Bilaplacian tests")
dx=dy=Interval()
d=dx*dy
Dx=Derivative(dx);Dy=Derivative(dy)
L=Dx^4⊗I+2*Dx^2⊗Dy^2+I⊗Dy^4


A=[dirichlet(dx)⊗eye(dy);
        eye(dx)⊗dirichlet(dy);
        neumann(dx)⊗eye(dy);
        eye(dx)⊗neumann(dy);
         L]


u=linsolve(A,[ones(4);zeros(5)];tolerance=1E-5)
@test_approx_eq u(0.1,0.2) 1.0


F=[Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[2]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[3]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[4]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[5]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[6]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[7]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[8]));
    0]

u=linsolve(A,F;tolerance=1E-10)

@test_approx_eq u(0.1,0.2)  exp(0.1)*cos(0.2)




A=[(ldirichlet(dx)+lneumann(dx))⊗eye(dy);
        (rdirichlet(dx)+rneumann(dx))⊗eye(dy);
        eye(dx)⊗(ldirichlet(dy)+lneumann(dy));
        eye(dx)⊗(rdirichlet(dy)+rneumann(dy));
        (ldirichlet(dx)-lneumann(dx))⊗eye(dy);
        (rdirichlet(dx)-rneumann(dx))⊗eye(dy);
        eye(dx)⊗(ldirichlet(dy)-lneumann(dy));
        eye(dx)⊗(rdirichlet(dy)-rneumann(dy));
         L]


u=linsolve(A,[ones(8);0];tolerance=1E-5)
@test_approx_eq u(0.1,0.2) 1.0



F=[2Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]));
    2Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[2]));
    Fun((x,y)->real(exp(x+1.0im*y))-imag(exp(x+1.0im*y)),rangespace(A[3]));
    Fun((x,y)->real(exp(x+1.0im*y))-imag(exp(x+1.0im*y)),rangespace(A[4]));
    0;
    0;
    Fun((x,y)->real(exp(x+1.0im*y))+imag(exp(x+1.0im*y)),rangespace(A[7]));
    Fun((x,y)->real(exp(x+1.0im*y))+imag(exp(x+1.0im*y)),rangespace(A[8]));
    0]

u=linsolve(A,F;tolerance=1E-10)

@test_approx_eq u(0.1,0.2)  exp(0.1)*cos(0.2)



## Test periodic x interval

println("    Periodic x Interval tests")


dθ=PeriodicInterval();dt=Interval(0,1.)
d=dθ*dt
ε=0.1
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
u0=Fun(θ->exp(-20θ^2),dθ,20)
@time u=linsolve([I⊗ldirichlet(dt);Dt-ε*Dθ^2-Dθ],[u0;0.];tolerance=1E-4)
@test_approx_eq_eps u(0.1,0.2) 0.3103472600253807 1E-2


# Transport equation
dθ=PeriodicInterval(-2.,2.);dt=Interval(0,1.)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
u0=Fun(θ->exp(-20θ^2),dθ)
A=Dt+Dθ

testbandedblockbandedoperator(A)

@time u=linsolve([I⊗ldirichlet(dt);Dt+Dθ],[u0;0.0];tolerance=1E-6)
@test_approx_eq_eps u(0.2,0.1) u0(0.1) 1E-6


## Small diffusoion


println("    Time evolution tests")

dx=Interval();dt=Interval(0,0.2)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x,t=Fun(dx*dt)


B=0.0
C=0.0
V=B+C*x
ε=0.1
f=Fun(x->exp(-30x^2),dx)
u=linsolve([timedirichlet(d);Dt-ε*Dx^2-V*Dx],[f;zeros(3)];tolerance=1E-6)

@test_approx_eq u(.1,.2) 0.496524222625512
B=0.1
C=0.2
V=B+C*x
u=linsolve([timedirichlet(d);Dt-ε*Dx^2-V*Dx],[f;zeros(3)];tolerance=1E-7)
@test_approx_eq u(.1,.2) 0.46810331039791464


## Periodic

println("    Periodic tests")

d=PeriodicInterval()^2
f=Fun((θ,ϕ)->exp(-10(sin(θ/2)^2+sin(ϕ/2)^2)),d)
A=Laplacian(d)+.1I
@time u=A\f
@test_approx_eq u(.1,.2) u(.2,.1)


d=PeriodicInterval()*Interval()
g=Fun((x,y)->real(cos(x+im*y)),∂(d))  # boundary data
@time u=[Dirichlet(d);Laplacian(d)]\Any[g;0.]

@test_approx_eq u(.1,.2) real(cos(.1+.2im))




dθ=PeriodicInterval(-2.,2.);dt=Interval(0,1.)


# Check bug in cache
CO=cache(ldirichlet(dt))
ApproxFun.resizedata!(CO,:,2)
ApproxFun.resizedata!(CO,:,4)
@test_approx_eq CO*Fun(exp,dt) 1.0


dθ=PeriodicInterval(-2.,2.);dt=Interval(0,3.)
d=dt*dθ
Dt=Derivative(d,[1,0]);Dθ=Derivative(d,[0,1])
A=[ldirichlet(dt)⊗I;Dt+Dθ]
testbandedblockbandedoperator(A[2])

u0=Fun(θ->exp(-20θ^2),dθ,20)
@time ut=linsolve(A,[u0;0.];tolerance=1E-5)
@test_approx_eq_eps ut(.1,.2) u0(.2-.1) 1E-6





# Beam


dθ=PeriodicInterval(0.0,1.0);dt=Interval(0,0.01)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1]);

B=[I⊗ldirichlet(dt),I⊗lneumann(dt)]
u0=Fun(θ->exp(-200(θ-.5).^2),dθ)
@time u=linsolve([B;Dt^2+Dθ^4],[u0;0.;0.];tolerance=1E-3)

@test_approx_eq_eps u(.1,.01) -0.2479768394633227  1E-3 #empirical

## Rectangle PDEs

println("    Rectangle tests")

# Screened Poisson

d=Interval()^2
@time u=linsolve([neumann(d);Laplacian(d)-100.0I],[ones(4);0.];tolerance=1E-12)
@test_approx_eq u(.1,.9) 0.03679861429138079

# PiecewisePDE

a=Fun((-1..1) \ [0,0.5],[1,0.5,1])
s=space(a)
dt=Interval(0,2.)
Dx=Derivative(s);Dt=Derivative(dt)
Bx=[ldirichlet(s);continuity(s,0)]

# test resize bug
CO=cache(Bx[2])
@test ApproxFun.colstop(CO.op,2) == 2
ApproxFun.resizedata!(CO,:,2)
ApproxFun.resizedata!(CO,:,4)
@test_approx_eq (CO*collect(1:4)).coefficients [3.,-1.]



## Test error


dx=Interval();dt=Interval(0,2.)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x,y=Fun(identity,d)
@time u=linsolve([I⊗ldirichlet(dt);Dt+x*Dx],[Fun(x->exp(-20x^2),dx);0.];tolerance=1E-12)

@test_approx_eq u(0.1,0.2) 0.8745340845783758  # empirical


dθ=PeriodicInterval();dt=Interval(0,1.)
d=dθ*dt
ε=0.1
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
u0=Fun(θ->exp(-20θ^2),dθ,20)
@time u=linsolve([I⊗ldirichlet(dt);Dt-ε*Dθ^2-Dθ],[u0;0.];tolerance=1E-4)
@test_approx_eq_eps u(0.1,0.2) 0.3103472600253807 1E-2
