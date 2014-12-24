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
println("Clenshaw large coeffs, 1 point: Time should be ~9e-6")

# @time is 8.853e-6 seconds


f=Fun(exp)
x=sample(f,100000)
x=sample(f,100000)
@time x=sample(f,100000)
println("Sample: Time should be ~0.27")
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds


f=Fun(x->cos(1000x),1000)
roots(f)
roots(f)
@time roots(f)
println("Roots: Time should be ~0.18")


## ODEs

d=Interval(-20000.,20000.)
x=Fun(identity,d)
u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
@time u=[dirichlet(d),diff(d)^2+I]\[1.,0.]
println("Cos/Sin: should be ~0.05")


d=Interval(-1000.,5.)
x=Fun(identity,d)
u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
@time u=[dirichlet(d),diff(d)^2-x]\[1.,0.]
println("Airy: should be ~0.05")

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
println("Poly: should be ~0.06")


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
println("Cos: should be ~0.08")


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
println("Laplace: should be ~0.04, 0.02")



d=Interval()^2
S=schurfact([neumann(d),lap(d)+100I],100)
u=S\ones(4)
u=S\ones(4)
@time u=S\ones(4)
println("Neumann Helmholtz: should be ~0.035")




d=Interval(-300.,5.)
x=Fun(identity,d)
A=Derivative(d)^2-x
u=null(A)
u=null(A)
@time u=null(A)
println("Null Airy: should be ~0.075")
