using ApproxFun, Base.Test, Base.Test
    import ApproxFun: testbandedblockbandedoperator, testblockbandedoperator

## Check operators
S=JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ=Laplacian(S)

testbandedblockbandedoperator(Δ)

u=Fun((x,y)->sin(π*x)*sin(π*y),S)
f=-2π^2*u


QR=qrfact(Δ)
ApproxFun.resizedata!(QR,:,1000)
@time v=QR\f
@test norm((u-v).coefficients)<100eps()


QR=qrfact(Δ)
ApproxFun.resizedata!(QR.R,:,100)
ApproxFun.resizedata!(QR.R,:,1000)
@time v=QR\f
@test norm((u-v).coefficients)<100eps()

QR=qrfact(Δ)
@time v=QR\f
@test norm((u-v).coefficients)<100eps()






## Rectangle PDE
dx=dy=Interval()
d=dx*dy
g=Fun((x,y)->exp(x)*cos(y),∂(d))

B=Dirichlet(d)
testblockbandedoperator(B)
testbandedblockbandedoperator(Laplacian(d)+0.0I)

A=[Dirichlet(d);Laplacian(d)+0.0I]

testblockbandedoperator(ApproxFun.interlace(A))

@time u=A\[g,0.]


@test u(.1,.2) ≈ real(exp(0.1+0.2im))


println("    Poisson tests")


## Poisson

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),Interval()^2,500)  #default is [-1,1]^2
d=domain(f)
A=[Dirichlet(d);Laplacian(d)]
@time  u=\(A,[zeros(∂(d));f];tolerance=1E-7)
@test_approx_eq_eps u(.1,.2) -0.04251891975068446 1E-5



# fourth order

println("    Bilaplacian tests")
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


@time u=\(A,[ones(4);zeros(5)];tolerance=1E-5)
@test u(0.1,0.2) ≈ 1.0




## Test periodic x interval

println("    Periodic x Interval tests")
d=PeriodicInterval()*Interval()

u_ex=Fun((x,y)->real(cos(x+im*y)),d)

B=Dirichlet(Space(d))

g=Fun((x,y)->real(cos(x+im*y)),rangespace(B))  # boundary data
@test norm((B*u_ex-g).coefficients) < 100eps()

testbandedblockbandedoperator(Laplacian(d))

@time u=[B;Laplacian(d)]\[g;0.]

@test u(.1,.2) ≈ real(cos(.1+.2im))



## Small diffusoion


println("    Time evolution tests")


## Schrodinger

dx=Interval(0.,1.);dt=Interval(0.0,0.001)
C=Conversion(Chebyshev(dx)*Ultraspherical(1,dt),Ultraspherical(2,dx)*Ultraspherical(1,dt))
testbandedblockbandedoperator(C)
testbandedblockbandedoperator(Operator{Complex128}(C))


d=dx*dt

x,y=Fun(d)
V=x^2

Dt=Derivative(d,[0,1]);Dx=Derivative(d,[1,0])

ϵ=1.
u0=Fun(x->exp(-100*(x-.5)^2)*exp(-1./(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
L=ϵ*Dt+(.5im*ϵ^2*Dx^2)
testbandedblockbandedoperator(L)

@time u=\([timedirichlet(d);L],[u0;zeros(3)];tolerance=1E-5)
@test u(0.5,0.001) 0.857215539785593+0.08694948835021317im  # empircal from ≈ schurfact


#
