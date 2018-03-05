using ApproxFun, Compat.Test
    import ApproxFun: Block, BlockBandedMatrix

## PDEs


S = JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ = Laplacian(S)

f = Fun((x,y)->sin(π*x)*sin(π*y),S)

QR=qrfact(Δ)
    ApproxFun.resizedata!(QR,:,400)
QR=qrfact(Δ)
    @time Δ[Block.(1:40), Block.(1:40)]
    @time ApproxFun.resizedata!(QR,:,400)
    @time \(QR,f; tolerance=1E-10)
println("Laplace Dirichlet: should be ~0.015, 0.015, 0.001")


d=Interval()^2
#dirichlet(d) is u[-1,:],u[1,:],u[:,-1],u[:,1]
A=[Dirichlet(d);Laplacian(d)]
f=Fun((x,y)->real(exp(x+im*y)),∂(d))


QR=qrfact(A)
    ApproxFun.resizedata!(QR,:,150)
    \(QR,[f;0.];tolerance=1E-10)


QR.R.data.block_sizes

hasmatchingblocks(A)

hasmatchingblocks(A)

QR.ncols
V = view(QR.R.data,1:153,1:153)


b = randn(153,153)
UpperTriangular(V) \ b

kr = jr = 1:153
N,  N_n = BlockBandedMatrices._find_block(BlockBandedMatrices.block_sizes(A), 1, kr[end])

view(V,
≈


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
