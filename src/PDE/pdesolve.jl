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



###
# pde_standardize_rhs manipulate "non-standard" rhs
# so it has a standardized form
###
function pde_standardize_rhs(A,f::Vector)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)

    rs=rangespace(A)

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
            F=zeros(rs)
        else
            F=Fun(f[end],rs)
        end
    else
        f=pad(f,max(length(indsBx)+length(indsBy),length(f)))
    
        fx=isempty(indsBx)?[]:convert2fun(f[indsBx],domainspace(A,2))
        fy=isempty(indsBy)?[]:convert2fun(f[indsBy],domainspace(A,1))       
        
        if length(f)==length(indsBx)+length(indsBy)+1
            F=Fun(f[end],rs)
        else
            F=zeros(eltype(A),rs)
        end
    end
    
    fx,fy,F
end   




function pde_standardize_rhs(A,f::Matrix)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)
    @assert length(indsBx)+length(indsBy)==size(f,1)


    fx=isempty(indsBx)?[]:f[indsBx,:]
    fy=isempty(indsBy)?[]:f[indsBy,:]
    F=fill(zeros(rangespace(A)),1,size(f,2))
    
    fx,fy,F
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


function pdesolve(A::AbstractPDEOperatorSchur,f::Array,nx=100000)
    fx,fy,F=pde_standardize_rhs(A,f)
    Fun(cont_constrained_lyap(A,fx,fy,F,nx),domainspace(A))
end

pdesolve(A::AbstractPDEOperatorSchur,f::Union(Fun,MultivariateFun,Number),nx...)=pdesolve(A,[f],nx...)
pdesolve{T<:PDEOperator}(A::Vector{T},f::Vector,n::Integer,n2...)=pdesolve(schurfact(A,n),f,n2...)
pdesolve{T<:PDEOperator}(A::Vector{T},f::Union(Fun,MultivariateFun,Number),n...)=pdesolve(A,[f],n...)
pdesolve(A::PDEOperator,f...)=pdesolve([A],f...)




\{T<:PDEOperator}(A::Vector{T},f::Array)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Array)=pdesolve(A,f)
\{T<:PDEOperator}(A::Vector{T},f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Fun)=pdesolve(A,f)
\(A::AbstractPDEOperatorSchur,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
\(A::PDEOperator,f::Union(MultivariateFun,Number,Fun))=pdesolve(A,f)
