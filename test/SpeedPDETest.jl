using ApproxFun, Compat.Test
    import ApproxFun: Block, BlockBandedMatrix

## PDEs


S = JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ = Laplacian(S)

f = Fun((x,y)->sin(π*x)*sin(π*y),S)

QR=qrfact(Δ)
    QR.R

BLAS.axpy!(1.0, view(B.op,KR,JR), view(B.data,KR,JR))
Δ[Block.(1:5), Block.(1:5)]
    ApproxFun.resizedata!(QR,:,400)
    \(QR,f;tolerance=1E-10)

@time Δ[Block.(1:40), Block.(1:40)]


V = view(Δ, Block.(1:40), Block.(1:40))
@time convert(BlockBandedMatrix, V)



Δ.op.ops[1]
@time Ao, Bo = (Δ.op.ops[1].ops)
@time A = Ao[Block.(1:40), Block.(1:40)];
@time B = Ao[Block.(1:40), Block.(1:40)];

@time A*B
@time Δ.op.ops[1][Block.(1:40), Block.(1:40)]
@time Δ.op[Block.(1:40), Block.(1:40)]





@time Δ.op.ops[1].ops[][Block.(1:40), Block.(1:40)]
@time Δ.op.ops[1].ops[3][Block.(1:40), Block.(1:40)]
@time Δ.op.ops[2][Block.(1:40), Block.(1:40)]


QR.R.datasize
QR=qrfact(Δ)
    @profile ApproxFun.resizedata!(QR,:,400)

Profile.print()

ProfileView.view()

using ProfileView


QR=qrfact(Δ)
    @time Δ[Block.(1:40), Block.(1:40)]
    @time ApproxFun.resizedata!(QR,:,400)
    @time \(QR,f; tolerance=1E-10)
using SO

@profile ApproxFun.resizedata!(QR,:,400)
Profile.clear()
@profile \(QR,f;tolerance=1E-10)

UpperTriangular(view(QR.R.data, 1:10, 1:10))

Profile.print()

println("Laplace Dirichlet: should be ~0.015, 0.001")


d=Interval()^2
#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
A=[Dirichlet(d);Laplacian(d)]
f=Fun((x,y)->real(exp(x+im*y)),∂(d))


QR=qrfact(A)
    ApproxFun.resizedata!(QR,:,150)
    \(QR,[f;0.];tolerance=1E-10)
QR=qrfact(A)
    @time ApproxFun.resizedata!(QR,:,150)
    @time u=\(QR,[f;0.];tolerance=1E-10)

println("Laplace: should be ~0.06, 0.001")



d=Interval()^2
[Neumann(d);lap(d)+100I]\[[[1,1],[1,1]],0]
@time [Neumann(d);lap(d)+100I]\[[[1,1],[1,1]],0]
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
