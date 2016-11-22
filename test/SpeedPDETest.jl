using ApproxFun, Base.Test


## PDEs


S=JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ=Laplacian(S)

f = Fun((x,y)->sin(π*x)*sin(π*y),S)
QR=qrfact(Δ)
    @time ApproxFun.resizedata!(QR,:,400)
    @time linsolve(QR,f;tolerance=1E-10)
QR=qrfact(Δ)
    @time ApproxFun.resizedata!(QR,:,400)
    @time linsolve(QR,f;tolerance=1E-10)



println("Laplace Dirichlet: should be ~0.05, 0.003")


d=Interval()^2


#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
A=[Dirichlet(d);Laplacian(d)]
f=Fun((x,y)->real(exp(x+im*y)),∂(d))

QR=qrfact(A)
    @time ApproxFun.resizedata!(QR.R,:,400)
    @time ApproxFun.resizedata!(QR,:,200)


@time linsolve(QR,[f;0.];tolerance=1E-10)
QR=qrfact(A)
    @time ApproxFun.resizedata!(QR,:,150)
    @time u=linsolve(QR,[f;0.];tolerance=1E-10)

println("Laplace: should be ~0.06, 0.001")



d=Interval()^2
@time [neumann(d);lap(d)+100I]\[ones(4);0]
@time [neumann(d);lap(d)+100I]\[ones(4);0]
println("Neumann Helmholtz: should be ~0.032")



#
# dx=Interval(0.,1.);dt=Interval(0.0,0.54)
# d=dx*dt
#
# x,y=Fun(d)
# V=x^2
#
# Dt=Derivative(d,[0,1]);Dx=Derivative(d,[1,0])
#
# ϵ=0.1
#
# u0=Fun(x->exp(-25*(x-.5)^2)*exp(-1.0im/(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
# L=1im*ϵ*Dt+.5*ϵ^2*Dx^2-V⊗1
#
# PO=discretize([timedirichlet(d);L],50)
# @time PO=discretize([timedirichlet(d);L],50)
# u=PO\u0
# @time    u=PO\u0
#
# println("Schrodinger: should be ~0.013,0.015")
