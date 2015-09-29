using ApproxFun, Base.Test

# This avoids getting killed on Travis.cl
if OS_NAME == :Darwin
    gc_enable(false)
end

## PDEs

d=Interval()^2

x=Fun(identity,d[1]);y=Fun(identity,d[2])

#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
A=[dirichlet(d);lap(d)]
f=[real(exp(-1+1.im*y));
                        real(exp(1+1im*y));
                        real(exp(x-1im));
                        real(exp(x+1im))]
S=schurfact(A,100)
@time S=schurfact(A,100)
u=S\f
u=S\f
@time u=S\f;
println("Laplace: should be ~0.014, 0.01")



d=Interval()^2
S=schurfact([neumann(d);lap(d)+100I],100)
@time S=schurfact([neumann(d);lap(d)+100I],100)
u=S\ones(4)
u=S\ones(4)
@time u=S\ones(4)
println("Neumann Helmholtz: should be ~0.016, 0.016")




dx=Interval(0.,1.);dt=Interval(0.0,0.54)
d=dx*dt

V=Fun(x->x^2,dx)

Dt=diff(d,2);Dx=diff(d,1)

ϵ=0.1

u0=Fun(x->exp(-25*(x-.5)^2)*exp(-1.im/(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
L=1im*ϵ*Dt+.5*ϵ^2*Dx^2-V⊗1

PO=discretize([timedirichlet(d);L],50)
@time PO=discretize([timedirichlet(d);L],50)
u=PO\u0
@time    u=PO\u0

println("Schrodinger: should be ~0.013,0.015")


if OS_NAME == :Darwin
    gc_enable(true)
end
