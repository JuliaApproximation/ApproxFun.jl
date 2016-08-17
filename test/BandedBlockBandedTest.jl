l=u=1
λ=μ=1
N=M=10
cols=rows=1:N
data=ones(λ+μ+1,(l+u+1)*sum(cols))
A=ApproxFun.BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols)
v=ones(size(A,1))
M=full(A)
@test norm(A*v-M*v) ≤ 100eps()
