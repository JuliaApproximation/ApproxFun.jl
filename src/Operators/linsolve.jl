##commondomainspace

##TODO: is this too hacky?
domainspace(f::Fun)=space(f)

function commondomainspace(P::Vector)
    ret = AnySpace()
    
    for op in P
        d = domainspace(op)
        @assert ret == AnySpace() || d == AnySpace() || typeof(ret) == typeof(d)
        
        if d != AnySpace()
            ret = d
            
            if domain(ret) != AnyDomain()
                return ret
            end
        end
    end
    
    ret
end

commondomainspace{T<:Number}(P::Vector,g::Array{T})=commondomainspace(P)
commondomainspace(P::Vector,g)=commondomainspace([P,g])



## Linear Solve




Fun_coefficients(b::Vector,sp)=vcat(map(f-> isa(f,Fun)? coefficients(f,sp) :  f,b)...)


function Fun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    A=promotedomainspace(A)
    #TODO: Use information from b to promote range space if necessary
    u=adaptiveqr(A,Fun_coefficients(b,rangespace(A[end])),tolerance,maxlength)  ##TODO: depends on ordering of A
    
    Fun(u,commondomainspace(A,b))
end

function Fun_linsolve{T<:Operator,N<:Number}(A::Vector{T},b::Array{N,2};tolerance=0.01eps(),maxlength=1000000)
    u=adaptiveqr(A,b,tolerance,maxlength)  ##TODO: depends on ordering of A
    d=commondomain(A)
    Fun[Fun(u[:,k],d) for k=1:size(u,2)]
end



# function FFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
#     @assert length(A) == 1
# 
#     u=adaptiveqr([interlace(A[1])],FFun_coefficients(b),tolerance,maxlength)
#     
#     FFun(deinterlace(u),commondomain(A,b))    
# end

function linsolve{T<:Operator}(A::Vector{T},b::Array;tolerance=0.01eps(),maxlength=1000000)
    d=commondomain(A,b)

    if typeof(d) <: Domain
        Fun_linsolve(A,b;tolerance=tolerance,maxlength=maxlength)
    else
        adaptiveqr(A,b,tolerance,maxlength)
    end    
end


 function linsolve{T<:Operator,M<:Number}(A::Array{T,2},b::Vector{M};tolerance=0.01eps(),maxlength=1000000)
    m = size(A,2)
 
     ret=adaptiveqr(interlace(A),b,tolerance,maxlength)  #Given just an array, we don't know how to interlace
                                                         #so assume user knows, this is correct for bc rows
                                     
                                     
     Fun[Fun(ret[k:m:end],commondomain(A[:,k])) for k=1:m]
 end
 
 
function linsolve{T<:Operator,M<:Number}(A::Array{T,2},b::Array{M,2};tolerance=0.01eps(),maxlength=1000000)
    m = size(A,2)
    
    ret=adaptiveqr(interlace(A),b,tolerance,maxlength)  #Given just an array, we don't know how to interlace
                                                     #so assume user knows, this is correct for bc rows
                                 
                                 
    [Fun(ret[k:m:end,j],commondomain(A[:,k])) for k=1:m,j=1:size(b,2)]
end
 

scalarorfuntype{S,T<:Number}(::Fun{S,T})=T
scalarorfuntype{T<:Number}(::T)=T
scalarorfuntype{T<:Number}(b::Vector{T})=T
scalarorfuntype(b::Vector{Any})=promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::Vector{F})=promote_type(map(scalarorfuntype,b)...)
 
function linsolve{T<:Operator}(A::Array{T,2},b::Vector{Any};kwds...)
    m,n=size(A)

    br=m-n

    l=mapreduce(length,max,b[br+1:end])

    r=zeros(scalarorfuntype(b),br+n*l)
    
    r[1:br]=b[1:br]
    
    for k=br+1:m
        sp=findmaxrangespace([A[k,:]...])
        if k > length(b)## assume its zer
            r[k:n:end]=zeros(l)
        elseif isa(b[k],Fun)
            ##TODO: boiunds check
            r[k:n:end]=pad(coefficients(b[k],sp),l)
        else  #type is scalar
            r[k:n:end]=pad([b[k]],l)
        end
    end

    linsolve(A,r;kwds...)
end 


linsolve(A::Operator,b::Array;kwds...)=linsolve([A],b;kwds...)
linsolve(A,b;kwds...)=linsolve(A,[b];kwds...)


\{T<:Operator}(A::Array{T,2},b::Array)=linsolve(A,b)
\{T<:Operator}(A::Vector{T},b::Array)=linsolve(A,b)
\(A::Operator,b)=linsolve(A,b)

