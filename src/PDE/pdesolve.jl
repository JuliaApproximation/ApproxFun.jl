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

function convert2funvec{T<:Number}(f::Vector{T},d::IntervalDomainSpace)
    ret=Array(Fun{T},length(f))
    for k=1:length(f)
        ret[k]=Fun(f[k],d)
    end
    ret
end
convert2funvec{T<:Fun}(f::Vector{T},d::IntervalDomainSpace)=f
function convert2funvec(f::Vector{Any},d::IntervalDomainSpace)
    mytyp=Fun{Float64}
    
    for fk in f
        if typeof(fk) == Fun{Complex{Float64}}
            mytyp=Fun{Complex{Float64}}
        end
    end
    
    ret=Array(mytyp,length(f))
    
    for k=1:length(f)
        ##TODO: Check domains match for Funs
        ret[k]=Fun(f[k],d)
    end
    ret
end


function pdesolve_mat(A::PDEOperatorSchur,f::Vector,nx=100000)
    if length(f) < length(A.indsBx)+length(A.indsBy)+1
        f=[f,zeros(length(A.indsBx)+length(A.indsBy)+1-length(f))]
    end

    ##TODO: makes more sense as a domain space of the boundary ops once thats set up
    fx=convert2funvec(f[A.indsBx],ChebyshevSpace(domain(A,2)))
    fy=convert2funvec(f[A.indsBy],ChebyshevSpace(domain(A,1)))
    

    ff=f[end]
    if typeof(ff)<:Number
        F=zeros(1,size(A.S,1)-numbcs(A.S)) 
        F[1,1]=ff
    elseif typeof(ff)<:Fun && domain(ff) == AnyDomain()
        ##TODO: beter method of telling constant fun
        F=zeros(1,size(A.S,1)-numbcs(A.S)) 
        F[1,1]=ff.coefficients[1]        
    else # typeof(ff) <:Fun2D || TensorFun
        F=coefficients(ff,rangespace(A,1),rangespace(A,2))
    end        
    

    cont_constrained_lyap(A,fy,fx,F,nx)
end

pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,nx::Integer,ny::Integer)=pdesolve_mat(PDEOperatorSchur(A,ny),f,nx)
pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,ny::Integer)=pdesolve_mat(PDEOperatorSchur(A,ny),f)



pdesolve_mat{T<:PDEOperator}(A::Vector{T},f)=pdesolve_mat(A,f,10000eps())
function pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,tol::Real)
    @assert tol>0
    maxit=11
   
    for k=5:maxit
        u=pdesolve_mat(A,f,2^k)
        if norm(map(f->norm(f.coefficients),u[end-2:end]))<tol
            return u
        end
    end
    error("Maximum number of iterations " * string(maxit) * "reached")
end



pdesolve(A::PDEOperatorSchur,f::Vector,nx...)=TensorFun(pdesolve_mat(A,f,nx...),domain(A,2))
pdesolve{T<:PDEOperator}(A::Vector{T},f::Vector)=TensorFun(pdesolve_mat(A,f),domain(A[end],2))
pdesolve{T<:PDEOperator}(A::Vector{T},f::Vector,n...)=TensorFun(pdesolve_mat(A,f,n...),domain(A[end],2))
pdesolve{T<:PDEOperator}(A::Vector{T},f::Fun,n...)=pdesolve(A,[f],n...)




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
\{T<:PDEOperator}(A::Vector{T},f::Fun)=pdesolve(A,f)
\(A::PDEOperatorSchur,f::Fun)=pdesolve(A,f)

