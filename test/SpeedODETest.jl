using ApproxFun, SpecialFunctions, Test
    import ApproxFun: ldiv_coefficients

## ODEs

d = Interval(-20000.,20000.)
x = Fun(identity,d)
u = [Dirichlet(d);Derivative(d)^2+I]\[[1,0],0]
u = [Dirichlet(d);Derivative(d)^2+I]\[[1,0],0]
@time u = [Dirichlet(d);Derivative(d)^2+I]\[[1,0],0]

println("Cos/Sin: should be ~0.016920 seconds (3.19 k allocations: 12.593 MiB)")

d=Interval(-1000.,5.)
x=Fun(identity,d)
u=[Dirichlet(d);Derivative(d)^2-x]\[[1,0],0]
u=[Dirichlet(d);Derivative(d)^2-x]\[[1,0],0]
@time u=[Dirichlet(d);Derivative(d)^2-x]\[[1,0],0]
println("Airy: 0.014356 seconds (1.08 k allocations: 8.015 MB)")

M=cache([Dirichlet(d);Derivative(d)^2-x];padding=true)
@time ApproxFun.resizedata!(M,12500,:)
println("Airy construct op: 0.003200 seconds (131 allocations: 3.723 MiB)")



S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+(7+2x+6x^2)
B=Dirichlet(S)
n=20000
rhs=fill(1.0,n+2)
u=ldiv_coefficients([B;L],rhs)
u=ldiv_coefficients([B;L],rhs)
@time u=ldiv_coefficients([B;L],rhs)
println("Poly: should be ~0.020926 seconds (3.00 k allocations: 9.416 MiB)")


S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+cos(x)
B=Dirichlet(S)
n=2000
rhs=fill(1.0,n+2)
u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
@time u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
println("Cos: should be ~0.0075")

S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+sin(x)
B=Dirichlet(S)
n=2000
rhs=fill(1.0,n+2)
u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
@time u=ldiv_coefficients([B;L],rhs;maxlength=Inf)
println("Sin: should be ~0.008663 seconds (660 allocations: 2.987 MB)")

## Bessel

x=Fun(identity,1..2000)
d=domain(x)
B=Dirichlet()
ν=1000.0
L=x^2*𝒟^2 + x*𝒟 + (x^2 - ν^2)   # our differential operator
u=[B;L]\[[besselj(ν,first(d)),besselj(ν,last(d))],0]
@time u=[B;L]\[[besselj(ν,first(d)),besselj(ν,last(d))],0]
println("Bessel: should be ~0.008441 seconds (6.14 k allocations: 4.765 MiB)")


x=Fun()
exp(10000*im*x)
@time exp(10000*im*x)
println("Complex exp: Time should be 0.03")


## Piecewise
x=Fun(identity,Domain(-20..15) \ Set([-10.,-5.,0.,1.]))
sp=space(x)
D=Derivative(sp)
B=[Dirichlet(sp);continuity(sp,0:1)]
u=[B;
    D^2-x]\[[airyai(-20.),0.],zeros(8),0];
@time u=[B;
    D^2-x]\[[airyai(-20.),0.],zeros(8),0]


println("Piecewise Airy: should be ~0.008")


## Vector
d=ChebyshevInterval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();

f=Fun(x->[exp(x),cos(x)],d)

A=[B 0;
   Bn 0;
   0 B;
   D^2-I 2.0I;
   0 D+I];

u=A\Any[0.,0.,0.,f...]
@time u=A\Any[0.,0.,0.,f...]
println("Systems: should be ~0.0008")


d=Interval(-300.,5.)
x=Fun(identity,d)
A=Derivative(d)^2-x
u=nullspace(A)
@test A[1:10,1:10] ≈ transpose(transpose(A)[1:10,1:10])
@time u=nullspace(A)
println("Nullspace Airy: 0.052730 seconds (75.21 k allocations: 56.736 MB)")
