export pdesolve


include("cont_lyap.jl")



# Converts a vector of constants/funs to a vector of Funs
function convert2fun{T<:Number,S<:Space}(f::Array{T},sp::S)
    ret=similar(f,Fun{S,T})
    for k=1:size(f,1),j=1:size(f,2)
        ret[k,j]=Fun(f[k,j],sp)
    end
    ret
end
convert2fun{T<:Fun}(f::Array{T},d::Space)=f
function convert2fun{S<:Space}(f::Array,sp::S)
    mytyp=Fun{S,promote_type(mapreduce(eltype,promote_type,f),eltype(sp))}

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


function depiecereorient(f,d)
    vf=pieces(f)
    dd1=d[1];dd2=d[2]

    if length(vf)==4
        # assume we are in boudnary
        # boundary has positive orientation
        vf=Any[setdomain(reverseorientation(vf[4]),dd2),
               setdomain(vf[2],dd2),
               setdomain(vf[1],dd1),
               setdomain(reverseorientation(vf[3]),dd1)]
    else
        @assert length(vf)==2
        if isa(dd1,PeriodicDomain)
            vf=Any[setdomain(vf[1],dd1),setdomain(reverseorientation(vf[2]),dd1)]
        elseif isa(dd2,PeriodicDomain)
            vf=Any[setdomain(vf[1],dd2),setdomain(reverseorientation(vf[2]),dd2)]
        end
    end
    vf
end

function pde_standardize_rhs(A,f::Vector)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)

    rs=rangespace(A)

    # if the fun lives on the boundary, we need to get rid of
    # the boundary information
    # i.e., if it lives on (-1-im,-1+im) we want to convert it to
    # living on (-1,1)
    if isa(A,PDEOperatorSchur) &&
            !isempty(f) &&
            isa(f[1],Fun) &&
            domain(f[1])==∂(domain(A))
        @assert length(f)<=2
        @assert isa(f,Vector)

        #TODO: More elegent domain conversion
        vf=depiecereorient(f[1],domain(A))

        if length(f)==1
            f=vf
        else
            @assert length(f)==2
            f=Any[vf...;f[end]]
        end
    end
    f=pad(f,max(length(indsBx)+length(indsBy),length(f)))

    ds=domainspace(A)
    fx=isempty(indsBx)?Fun{spacetype(ds,2),eltype(A)}[]:convert2fun(f[indsBx],ds[2])
    fy=isempty(indsBy)?Fun{spacetype(ds,1),eltype(A)}[]:convert2fun(f[indsBy],ds[1])

    if length(f)==length(indsBx)+length(indsBy)+1
        F=Fun(f[end],rs)
    else
        F=zeros(eltype(A),rs)
    end

    fx,fy,F
end




function pde_standardize_rhs(A,f::Matrix)
    ∂A=∂(domain(A))
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)

    if isa(A,PDEOperatorSchur) &&
                !isempty(f) &&
                all(g->isa(g,Fun)&&domain(g)==∂A,f[1,:])
        @assert size(f,1)==1
        d=domain(A)
        f2=Array(Any,length(indsBx)+length(indsBy),size(f,2))
        for k=1:size(f,2)
            f2[:,k]=depiecereorient(f[1,k],d)
        end
        f=f2
    end

    @assert length(indsBx)+length(indsBy)==size(f,1)


    fx=isempty(indsBx)?[]:convert2fun(f[indsBx,:],domainspace(A,2))
    fy=isempty(indsBy)?[]:convert2fun(f[indsBy,:],domainspace(A,1))
    F=fill(zeros(rangespace(A)),1,size(f,2))

    fx,fy,F
end


function pdesolve{T<:Operator}(A::Vector{T},f;tolerance::Float64=100000eps(),opts...)
    @assert tolerance>0
    maxit=11

    for k=5:maxit
        u=pdesolve(A,f,2^k;tolerance=tolerance/100,opts...)
        if norm(map(f->norm(f.coefficients),u.coefficients[end-2:end]))<tolerance
            return u
        end
    end
    error("Maximum number of iterations " * string(maxit) * "reached")
end

function pdesolve(A::AbstractPDEOperatorSchur,f::Vector;opts...)
    fx,fy,F=pde_standardize_rhs(A,f)
    ProductFun(cont_constrained_lyap(A,fx,fy,F;opts...),domainspace(A))
end


function pdesolve(A::AbstractPDEOperatorSchur,f::Matrix;opts...)
    fx,fy,F=pde_standardize_rhs(A,f)
    X=cont_constrained_lyap(A,fx,fy,F;opts...)
    ds=domainspace(A)
    #TODO: Change Float64
    ProductFun{typeof(ds[1]),
               typeof(ds[2]),
               typeof(ds),
               Float64}[ProductFun(X[:,k],ds) for k=1:size(X,2)]
end

pdesolve(A::AbstractPDEOperatorSchur,f::Union{Fun,MultivariateFun,Number},nx...;opts...)=pdesolve(A,[f],nx...;opts...)
pdesolve{T<:Operator}(A::Vector{T},f::Vector,n::Integer,n2...;opts...)=pdesolve(schurfact(A,n),f,n2...;opts...)
pdesolve{T<:Operator}(A::Vector{T},f::Union{Fun,MultivariateFun,Number},n...;opts...)=pdesolve(A,[f],n...;opts...)
pdesolve(A::Operator,f...;opts...)=pdesolve([A],f...;opts...)


\(A::AbstractPDEOperatorSchur,f::Union{MultivariateFun,Number,Fun,Array})=pdesolve(A,f)
