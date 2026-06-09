using ApproxFun, Base.Test


#BandedBlockBandedMatrix
l=u=1
λ=μ=1
N=M=10
cols=rows=1:N
data=ones(λ+μ+1,(l+u+1)*sum(cols))
A=ApproxFun.BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols)
v=ones(size(A,1))
M=full(A)
@test norm(A*v-M*v) ≤ 100eps()


A=ApproxFun.bbbrand(Float64,1,1,1,1,1:10,1:10)
B=ApproxFun.bbbrand(Float64,1,1,1,1,1:10,1:10)
@test_approx_eq A*B full(A)*full(B)


#Ragged Matrix

cols=Int[rand(1:k+2) for k=1:5]
B=ApproxFun.rrand(Float64,maximum(cols),cols)
cols=Int[rand(1:k+2) for k=1:size(B,1)]
A=ApproxFun.rrand(Float64,maximum(cols),cols)
@test_approx_eq full(A)*full(B) full(A*B)


## BandedBlockMatrix

N=10
A=ApproxFun.bbrand(Float64,1,1,1:N,1:N)
B=ApproxFun.bbrand(Float64,1,1,1:N,1:N)

@test_approx_eq full(A)*full(B) full(A*B)

v=ones(size(A,1))

M=full(A)
@test norm(A*v-M*v) ≤ 100eps()
