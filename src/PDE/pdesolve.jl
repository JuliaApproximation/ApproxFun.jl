export pdesolve


include("cont_lyap.jl")




# function pdetoarray(Bxin,Byin,Lin,Min,nx::Integer,ny::Integer)
#     Bx,Gx,Lx,Mx=pdetoarray(Bxin,Lin[1],Min[1],nx)
#     By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny)    
#     
# 
#     Gx=toarray(Gx,ny);Gy=toarray(Gy,nx)    
#     
#     
#     Bx,Gx,By,Gy,Lx,Ly,Mx,My    
# end

function convert2funvec{T<:Number,D<:IntervalDomain}(f::Vector{T},d::D)
    ret=Array(IFun{T,D},length(f))
    for k=1:length(f)
        ret[k]=IFun(f[k],d)
    end
    ret
end
convert2funvec{T<:IFun}(f::Vector{T},d::IntervalDomain)=f
function convert2funvec{D<:IntervalDomain}(f::Vector{Any},d::D)
    mytyp=IFun{Float64,D}
    
    for fk in f
        if typeof(fk) == IFun{Complex{Float64},D}
            mytyp=IFun{Complex{Float64},D}
        end
    end
    
    ret=Array(mytyp,length(f))
    
    for k=1:length(f)
        ##TODO: Check domains match for IFuns
        ret[k]=IFun(f[k],d)
    end
    ret
end


function pdesolve_mat(A::PDEOperatorSchur,f::Vector)
    if length(f) < length(A.indsBx)+length(A.indsBy)+1
        f=[f,zeros(length(A.indsBx)+length(A.indsBy)+1-length(f))]
    end

    fx=convert2funvec(f[A.indsBx],domain(A,2))
    fy=convert2funvec(f[A.indsBy],domain(A,1))
    

    ff=f[end]
    if typeof(ff)<:Number
        F=zeros(1,size(A.S,1)-numbcs(A.S)) 
        F[1,1]=ff
    else # typeof(ff) <:Fun2D || TensorFun
        F=coefficients(ff,rangespace(A,1).order,rangespace(A,2).order)
    end        
    

    cont_constrained_lyap(A,fy,fx,F)
end

pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,ny::Integer)=pdesolve_mat(PDEOperatorSchur(A,ny),f)



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



pdesolve(A::PDEOperatorSchur,f::Vector)=TensorFun(pdesolve_mat(A,f),domain(A,2))
pdesolve{T<:PDEOperator}(A::Vector{T},f)=TensorFun(pdesolve_mat(A,f),domain(A[end],2))
pdesolve{T<:PDEOperator}(A::Vector{T},f,ny)=TensorFun(pdesolve_mat(A,f,ny),domain(A[end],2))




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
\(A::PDEOperatorSchur,f::Vector)=pdesolve(A,f)
\{T<:PDEOperator}(A::Vector{T},f::IFun)=pdesolve(A,[f])
\(A::PDEOperatorSchur,f::IFun)=pdesolve(A,[f])

