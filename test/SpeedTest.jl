using ApproxFun, Base.Test

gc_disable()

c = rand(1000)
x=rand(10000)
f=Fun(c)
y=f[x]
y=f[x]

@time y=f[x]
println("Clenshaw large coeffs, many points: Time should be ~0.024")
# 0.012482274  with unsafe_view
# 0.024306262 with inbounds

y=f[.1]
y=f[.1]
y=f[.1]

@time y=f[.1];
println("Clenshaw large coeffs, 1 point: Time should be ~6e-6")

# @time is 8.853e-6 seconds


f=Fun(exp)
x=sample(f,100000)
x=sample(f,100000)
@time x=sample(f,100000)
println("Sample: Time should be ~0.25")
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds


f=Fun(x->cos(1000x),1000)
roots(f)
roots(f)
@time roots(f)
println("Roots: Time should be ~0.15")


## ODEs

d=Interval(-20000.,20000.)
x=Fun(identity,d)
u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
@time u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
println("Cos/Sin: should be ~0.017")


d=Interval(-1000.,5.)
x=Fun(identity,d)
u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
@time u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
println("Airy: should be ~0.013")

S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+(7+2x+6x^2)
B=dirichlet(S)
n=20000
rhs=ones(n+2)
u=[B,L]\rhs
u=[B,L]\rhs
@time u=[B,L]\rhs
println("Poly: should be ~0.025")


S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+cos(x)
B=dirichlet(S)
n=2000
rhs=ones(n+2)
u=linsolve([B,L],rhs;maxlength=Inf)
u=linsolve([B,L],rhs;maxlength=Inf)
@time u=linsolve([B,L],rhs;maxlength=Inf)
println("Cos: should be ~0.0075")

S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+sin(x)
B=dirichlet(S)
n=2000
rhs=ones(n+2)
u=linsolve([B,L],rhs;maxlength=Inf)
u=linsolve([B,L],rhs;maxlength=Inf)
@time u=linsolve([B,L],rhs;maxlength=Inf)
println("Sin: should be ~0.015")


## Piecewise

x=Fun(identity,[-20.,-10.,-5.,0.,1.,15.])
sp=space(x)
D=Derivative(sp)

u=[dirichlet(sp),
    D^2-x]\[airyai(-10.)];
@time u=[dirichlet(sp),
    D^2-x]\[airyai(-10.)];    

println("Piecewise Airy: should be ~0.016")


## Vector 
d=Interval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();

f=Fun(x->[exp(x),cos(x)],d)

A=[B 0;
   Bn 0;
   0 B;
   D^2-I 2.I;
   0 D+I];
   
u=A\Any[0.,0.,0.,f]
@time u=A\Any[0.,0.,0.,f]
println("Systems: should be ~0.0013")

## PDEs

d=Interval()^2

x=Fun(identity,d[1]);y=Fun(identity,d[2])

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
A=[dirichlet(d),lap(d)]
f=[real(exp(-1+1.im*y)),
                        real(exp(1+1im*y)),
                        real(exp(x-1im)),
                        real(exp(x+1im))]
u=pdesolve(A,f,100)
u=pdesolve(A,f,100)
@time u=pdesolve(A,f,100);
S=schurfact(A,100)
u=S\f
u=S\f
@time u=S\f;
println("Laplace: should be ~0.03, 0.011")



d=Interval()^2
S=schurfact([neumann(d),lap(d)+100I],100)
u=S\ones(4)
u=S\ones(4)
@time u=S\ones(4)
println("Neumann Helmholtz: should be ~0.017")




d=Interval(-300.,5.)
x=Fun(identity,d)
A=Derivative(d)^2-x
u=null(A)
u=null(A)
@time u=null(A)
println("Null Airy: should be ~0.061")
