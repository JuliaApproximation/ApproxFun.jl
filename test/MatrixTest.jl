using ApproxFun, Base.Test
    import ApproxFun: Block



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

# Tests bug in Complex
ret=ApproxFun.bbbzeros(Float64,0,4,0,4,[1,2,2],[1,2,2])
view(ret,Block(3),Block(3))[1,1]=2.0
@test ret[4,4]==2
@test ret[3,5]==0

ret=ApproxFun.bbbzeros(Float32,0,4,0,4,[1,2,2],[1,2,2])
view(ret,Block(3),Block(3))[1,1]=2.0
@test ret[4,4]==2
@test ret[3,5]==0



ret=ApproxFun.bbbzeros(Complex128,0,4,0,4,[1,2,2],[1,2,2])
view(ret,Block(3),Block(3))[1,1]=2.0
@test ret[4,4]==2
@test ret[3,5]==0


#Ragged Matrix

cols=Int[rand(1:k+2) for k=1:5]
B=ApproxFun.rrand(Float64,maximum(cols),cols)
cols=Int[rand(1:k+2) for k=1:size(B,1)]
A=ApproxFun.rrand(Float64,maximum(cols),cols)
@test_approx_eq full(A)*full(B) full(A*B)


## BandedBlockMatrix
N=10
A=ApproxFun.bbones(Float64,1,1,1:N,1:N)
@test A[1,1] == 1
A[1,1]=2
@test A[1,1] == 2
A[10,10]=2
@test A[10,10] == 2


N=10
A=ApproxFun.bbrand(Float64,1,1,1:N,1:N)
B=ApproxFun.bbrand(Float64,1,1,1:N,1:N)

@test_approx_eq full(A)*full(B) full(A*B)

v=ones(size(A,1))

M=full(A)
@test norm(A*v-M*v) ≤ 100eps()

# Check bug
A=ApproxFun.bbzeros(Float64,4,1,[4],[1])
A[2,1] = 3.0



## view

N=10
A=ApproxFun.bbones(Float64,1,1,1:N,1:N)

view(A,Block(1),Block(1))[1,1]=2
@test A[1,1] == 2
view(A,Block(1)[1:1],Block(1))[1,1]=3
@test A[1,1] == 3
view(A,Block(1),Block(1)[1:1])[1,1]=4
@test A[1,1] == 4
view(A,Block(1)[1:1],Block(1)[1:1])[1,1]=5
@test A[1,1] == 5

strides(view(A,Block(2),Block(3)))  == (1,2)


@test view(view(A,Block(4),Block(5)),1:2,1:3) == ones(2,3)
@test view(A,Block(4)[2:3],Block(5)[3:5]) == ones(2,3)


view(A,Block(3),Block(3))[2,1] = 2

@test unsafe_load(pointer(view(A,Block(3)[2:3],Block(3)[1:3]))) == 2
view(A,Block(3),Block(3))[2,2] = 3
@test unsafe_load(pointer(view(A,Block(3)[2:3],Block(3)[2:3]))) == 3


view(A,Block(3),Block(4))[2,3]=6
@test view(A,Block(3)[2:3],Block(4)[3:4])[1] == 6
@test unsafe_load(pointer(view(A,Block(3)[2:3],Block(4)[3:4]))) == 6

N=10
A=ApproxFun.bbones(Float64,1,1,1:N,1:N)

@test view(view(A,Block(3),Block(4)),2:3,:) == ones(2,4)
@test view(view(A,Block(3),Block(4)),:,2:3) == ones(3,2)
