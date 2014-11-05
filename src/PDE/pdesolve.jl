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

function convert2funvec{T<:Number,D}(f::Vector{T},d::D)
    ret=Array(Fun{D,T},length(f))
    for k=1:length(f)
        ret[k]=Fun(f[k],d)
    end
    ret
end
convert2funvec{T<:Fun,D}(f::Vector{T},d::D)=f
function convert2funvec{D}(f::Vector,d::D)
    mytyp=Fun{D,mapreduce(eltype,promote_type,f)}
    
    ret=Array(mytyp,length(f))
    
    for k=1:length(f)
        #TODO: Check domains match for Funs
        ret[k]=Fun(f[k],d)
    end
    ret
end

function pde_normalize_rhs(A,f)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)

    # vec if boundary fun
    if domain(f[1])==∂(domain(A))
        @assert length(f)<=2

        #TODO: general spaces, will break later though if inconsistent
        vf=vec(f[1])
        for g in vf
            @assert isa(space(g),ChebyshevSpace) 
        end
        
        fx=map(g->Fun(g.coefficients,domain(A,2)),vf[indsBx])
        fy=map(g->Fun(g.coefficients,domain(A,1)),vf[indsBy])   
        
        if length(f)<2
            ff=0.0
        else
            ff=f[end]
        end
    else
        ##TODO: makes more sense as a domain space of the boundary ops once thats set up    
        fx=isempty(indsBx)?[]:convert2funvec(f[indsBx],domainspace(A,2))
        fy=isempty(indsBy)?[]:convert2funvec(f[indsBy],domainspace(A,1))       
        
        if length(f) < length(indsBx)+length(indsBy)+1
            ff=0.0
        else
            ff=f[end]
        end      
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


function pdesolve_mat(A::AbstractPDEOperatorSchur,f::Vector,nx=100000)
    fx,fy,F=pde_normalize_rhs(A,f)
    #TODO: swap fy and fx order below
    cont_constrained_lyap(A,fy,fx,F,nx)
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


pdesolve(A::AbstractPDEOperatorSchur,f::Vector,nx...)=Fun(pdesolve_mat(A,f,nx...),domainspace(A))
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



\{T<:PDEOperator}(A::Vector{T},f::Vector)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Vector)=pdesolve(A,f)
\{T<:PDEOperator}(A::Vector{T},f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Fun)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::PDEOperator,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
