export pdesolve


include("cont_lyap.jl")



# Converts a vector of constants/funs to a vector of Funs
function convert2fun{T<:Number,S<:FunctionSpace}(f::Array{T},sp::S)
    ret=similar(f,Fun{S,T})
    for k=1:size(f,1),j=1:size(f,2)
        ret[k,j]=Fun(f[k,j],sp)
    end
    ret
end
convert2fun{T<:Fun}(f::Array{T},d::FunctionSpace)=f
function convert2fun{S<:FunctionSpace}(f::Array,sp::S)
    mytyp=Fun{S,mapreduce(eltype,promote_type,f)}
    
    ret=similar(f,mytyp)
    
    for k=1:size(f,1),j=1:size(f,2)
        #TODO: Check domains match for Funs
        ret[k,j]=Fun(f[k,j],sp)
    end
    ret
end

function pde_normalize_rhs(A,f::Vector)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)

    # if the fun lives on the boundary, we need to get rid of
    # the boundary information
    # i.e., if it lives on (-1-im,-1+im) we want to convert it to 
    # living on (-1,1)
    if isa(A,PDEOperatorSchur) && !isempty(f) && isa(f[1],Fun) && domain(f[1])==∂(domain(A))
        @assert length(f)<=2
        @assert isa(f,Vector)

        #TODO: More elegent domain conversion
        vf=pieces(f[1])
        ds1=domainspace(A,1);ds2=domainspace(A,2)
        
        fx=map(g->(@assert isa(space(g),typeof(ds2));
        Fun(g.coefficients,ds2)),vf[indsBx])
        fy=map(g->(@assert isa(space(g),typeof(ds1));
        Fun(g.coefficients,ds1)),vf[indsBy])   
        
        if length(f)<2
            ff=0.0
        else
            ff=f[end]
        end
    else
        f=pad(f,length(indsBx)+length(indsBy)+1)
        
        fx=isempty(indsBx)?[]:convert2fun(f[indsBx],domainspace(A,2))
        fy=isempty(indsBy)?[]:convert2fun(f[indsBy],domainspace(A,1))       
        
        ff=f[end]
    end
    

    if isa(ff,Number)
        F=zeros(typeof(ff),1,1) 
        F[1,1]=ff
    elseif isa(ff,Fun) && domain(ff) == AnyDomain()
        ##TODO: beter method of telling constant fun
        F=zeros(typeof(ff.coefficients[1]),1,1) 
        F[1,1]=ff.coefficients[1]        
    else # typeof(ff) <:LowRankFun || TensorFun
        F=coefficients(ff,rangespace(A))
    end     
    
    fx,fy,F
end   


function pdesolve_mat(A::AbstractPDEOperatorSchur,f::Array,nx=100000)
    fx,fy,F=pde_normalize_rhs(A,f)
    cont_constrained_lyap(A,fx,fy,F,nx)
end

pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,nx::Integer,ny::Integer)=pdesolve_mat(schurfact(A,ny),f,nx)
pdesolve_mat{T<:PDEOperator}(A::Vector{T},f,ny::Integer)=pdesolve_mat(schurfact(A,ny),f)


##TODO: Do we need pdesolve_mat adaptive routine?
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



pdesolve{T<:PDEOperator}(A::Vector{T},f)=pdesolve(A,f,10000eps())
function pdesolve{T<:PDEOperator}(A::Vector{T},f,tol::Real)
    @assert tol>0
    maxit=11
   
    for k=5:maxit
        u=pdesolve(A,f,2^k)
        if norm(map(f->norm(f.coefficients),u.coefficients[end-2:end]))<tol
            return u
        end
    end
    error("Maximum number of iterations " * string(maxit) * "reached")
end


function pdesolve(S::PDEStrideOperatorSchur,f::Vector,nx=100000)
    @assert length(f)≤5

    f=pad(f,5)
    uo=pdesolve(S.odd,[f[1],f[2],f[3]+f[4],f[end]],nx) 
    ue=pdesolve(S.even,[f[1],f[2],f[4]-f[3],f[end]],nx)     
    
    ret=Array(typeof(first(uo.coefficients)),length(uo.coefficients)+length(ue.coefficients))
    ret[1:2:end]=uo.coefficients
    ret[2:2:end]=ue.coefficients
    
    TensorFun(ret,domainspace(S.odd,2).space)
end

pdesolve(A::AbstractPDEOperatorSchur,f::Array,nx...)=Fun(pdesolve_mat(A,f,nx...),domainspace(A))
pdesolve(A::AbstractPDEOperatorSchur,f::Union(Fun,MultivariateFun,Number),nx...)=pdesolve(A,[f],nx...)
pdesolve{T<:PDEOperator}(A::Vector{T},f::Vector,n::Integer,n2...)=pdesolve(schurfact(A,n),f,n2...)
pdesolve{T<:PDEOperator}(A::Vector{T},f::Union(Fun,MultivariateFun,Number),n...)=pdesolve(A,[f],n...)
pdesolve(A::PDEOperator,f...)=pdesolve([A],f...)



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



\{T<:PDEOperator}(A::Vector{T},f::Array)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Array)=pdesolve(A,f)
\{T<:PDEOperator}(A::Vector{T},f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Fun)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::PDEOperator,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
