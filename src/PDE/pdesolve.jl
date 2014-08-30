export pdesolve


include("cont_lyap.jl")


##TODO: Unify with coefficients()
toarray{T<:Functional}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
toarray{T<:IFun}(B::Array{T},n)=Float64[    j<=length(B[k])?B[k].coefficients[j]:0 for  k=1:length(B),j=1:n];
function toarray(B::Array,n)
    ret = zeros(length(B),n)
    
    for k=1:length(B), j=1:n
        if  typeof(B[k]) <: IFun
            ret[k,j] = j<=length(B[k])?B[k].coefficients[j]:0 
        elseif typeof(B[k]) <: Number && j == 1
            ret[k,j] = B[k]
        end
    end

    ret
end

function toarray{T<:Operator}(A::Vector{T},n::Integer,m::Integer)
    ret = zeros(n,m)
    
    nbc = typeof(A[end])<:Functional?length(A):length(A)-1
    for k=1:nbc
        ret[k,:]=A[k][1:m]
    end
    
    if nbc < length(A)
        ret[nbc+1:end,:]=A[end][1:n-nbc,1:m]
    end
    
    ret
end

function pdetoarray(Byin,Lin,Min,ny::Integer)
    Yop=promotespaces([Lin,Min])

    
    By=toarray(Byin[1],ny)   
    nbcy=length(Byin[1])
    
    Ly=Yop[1][1:ny-nbcy,1:ny]
    My=Yop[2][1:ny-nbcy,1:ny]    
    
    By,Byin[2],Ly,My    
end

function pdetoarray(Bxin,Byin,Lin,Min,nx::Integer,ny::Integer)
    Bx,Gx,Lx,Mx=pdetoarray(Bxin,Lin[1],Min[1],nx)
    By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny)    
    

    Gx=toarray(Gx,ny);Gy=toarray(Gy,nx)    
    
    
    Bx,Gx,By,Gy,Lx,Ly,Mx,My    
end


isxfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,1])<:Functional
isyfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,2])<:Functional
ispdeop(B::PDEOperator)=!isxfunctional(B)&&!isyfunctional(B)


function pdesolve_mat{T<:PDEOperator}(A::Vector{T},f::Array,ny::Integer)
    inds=find(isxfunctional,A)
    Bx=Functional[(@assert Ai.ops[1,2]==ConstantOperator{Float64}(1.0); Ai.ops[1,1]) for Ai in A[inds]]
    fx=f[inds]
    inds=find(isyfunctional,A)
    By=Functional[(@assert Ai.ops[1,1]==ConstantOperator{Float64}(1.0); Ai.ops[1,2]) for Ai in A[inds]]


    fy=convert(Vector{IFun{Float64, Interval{Float64}}},f[inds])
    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)
    LL=A[end]
    ff=f[end]
    
    if typeof(ff)<:Number
        F=zeros(1,ny-length(By)) 
    else
        error("General RHS not implemented")
    end

    if size(LL.ops)==(2,2)
        cont_constrained_lyap({Bx,fx},{By,fy},
                                    {LL.ops[1,1],LL.ops[1,2]},
                                    {LL.ops[2,1],LL.ops[2,2]},F,ny)
    else
        error("Higher rank PDEs not implemented") 
    end
end



function pdesolve_mat{T<:PDEOperator}(A::Vector{T},f)
    maxit=11
    tol=10000eps()
    for k=5:maxit
        u=pdesolve_mat(A,f,2^k)
        if norm(map(f->norm(f.coefficients),u[end-2:end]))<tol
            return u
        end
    end
    error("Maximum number of iterations " * string(maxit) * "reached")
end

pdesolve{T<:PDEOperator}(A::Vector{T},f)=Fun2D(pdesolve_mat(A,f),domain(A[end],2))


# 
# function pdesolve(Bxin,Byin,Lin,Min,Fin::Number,nx::Integer,ny::Integer)
#     F=zeros(nx-length(Bxin[1]),ny-length(Byin[1]))
#     F[1,1]=Fin
#     
#     pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
# end
# 
# function pdesolve(Bxin,Byin,Lin,Min,Fin::Fun2D,nx::Integer,ny::Integer)
#     Xop=promotespaces([Lin[1],Min[1]])
#     Yop=promotespaces([Lin[2],Min[2]])
#     
#     Xsp=rangespace(Xop[1])
#     Ysp=rangespace(Yop[1])    
#     
#     nbcx=length(Bxin[1]);nbcy=length(Byin[1])    
#     
#     F=pad(coefficients(Fin,Xsp,Ysp),nx-nbcx,ny-nbcy)
#     pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
# end
# 
# 
# 
# function pdesolve(Bxin,Byin,Lin,Min,F::Array,nx::Integer,ny::Integer)
# 
#     Bx,Gx,By,Gy,Lx,Ly,Mx,My=pdetoarray(Bxin,Byin,Lin,Min,nx,ny)
#     
#     
#     Fun2D(constrained_lyap({Bx Gx; By Gy},{Lx,Ly},{Mx,My},F),
#             domain(Lin[1]),domain(Lin[2]))
# end
# 
# 
# function pdesolve(Bxin,Byin,Lin,Min,Fin::Number,ny::Integer)
#     F=zeros(1,ny-length(Byin[1]))
#     F[1,1]=Fin
#     
#     pdesolve(Bxin,Byin,Lin,Min,F,ny)
# end
# 
# pdesolve(Bxin,Byin,Lin,Min,Fin::Fun2D,ny::Integer)=pdesolve(Bxin,Byin,Lin,Min,coefficients(Fin,rangespace(Lin[1]+Min[1]).order,rangespace(Lin[2]+Min[2]).order),ny)
# 
# 
# 
# pdesolve(Bx,By,Lin,Min,F::Array,ny::Integer)=Fun2D(cont_constrained_lyap(Bx, By,Lin,Min,F,ny),domain(Lin[2]))



\{T<:PDEOperator}(A::Vector{T},f::Vector)=pdesolve(A,f)

