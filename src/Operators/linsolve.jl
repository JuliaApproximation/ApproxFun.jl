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


# function Fun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
# 
# end

# function Fun_linsolve{T<:Operator,N<:Number}(A::Vector{T},b::Array{N,2};tolerance=0.01eps(),maxlength=1000000)
#     u=adaptiveqr(A,b,tolerance,maxlength)  ##TODO: depends on ordering of A
#     d=commondomain(A)
#     Fun[Fun(u[:,k],d) for k=1:size(u,2)]
# end



# function FFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
#     @assert length(A) == 1
# 
#     u=adaptiveqr([interlace(A[1])],FFun_coefficients(b),tolerance,maxlength)
#     
#     FFun(deinterlace(u),commondomain(A,b))    
# end


function linsolve{T<:Operator,N<:Number}(A::Vector{T},b::Array{N};tolerance=0.01eps(),maxlength=1000000)
    Ad=promotedomainspace(A)
    ds=commondomainspace(Ad)
    r=adaptiveqr(Ad,b,tolerance,maxlength)
    isa(ds,AnySpace)?r:Fun(r,ds)
end

function linsolve{T<:Operator}(A::Vector{T},b::Array{Any};tolerance=0.01eps(),maxlength=1000000)
    @assert size(b,2)==1  #TODO: reimplemnt
    A=promotedomainspace(A)
    
    for k=1:length(A)-1
        @assert isa(A[k],Functional)
    end
    
    for k=1:min(length(A)-1,length(b))
        @assert isa(b[k],Number)
    end    
    
    #TODO: Use information from b to promote range space if necessary
    if length(b)<size(A,1)
        # the ... converts b to a tuple of numbers so that r is a number Vec    
        r=[b...]
    elseif length(b)==size(A,1)
        r=[b[1:end-1]...,coefficients(b[end],rangespace(A[end]))]
    else 
        # we have list of possible funs, devec
        #TODO: constants
        r=[b[1:size(A,1)-1]...,coefficients(devec(b[size(A,1):end]),rangespace(A[end]))]
    end
    
    u=adaptiveqr(A,r,tolerance,maxlength)  ##TODO: depends on ordering of A
    
    Fun(u,commondomainspace(A,b))
end

function linsolve{T<:Operator,F<:Fun}(A::Vector{T},b::Array{F};kwds...)
    r=Array(Any,size(b))
    
    # convert constant funs to constants
    # this undoes the effect of [0.,f]
    for k=1:size(A,1)-1,j=1:size(b,2)
        # we only allow constants
        @assert length(b[k,j])==1
        #TODO: 1,1 entry may not be zero
        r[k,j]=b[k,j].coefficients[1]
    end
    
    r[size(A,1):end,:]=b[size(A,1):end,:]
    
    linsolve(A,r;kwds...)
end

 
linsolve{T<:Operator}(A::Array{T,2},b;kwds...)=linsolve(interlace(A),b;kwds...)


linsolve(A::Operator,b::Array;kwds...)=linsolve([A],b;kwds...)
linsolve(A,b;kwds...)=linsolve(A,[b];kwds...)


\{T<:Operator}(A::Array{T,2},b::Array)=linsolve(A,b)
\{T<:Operator}(A::Vector{T},b::Array)=linsolve(A,b)
\(A::Operator,b)=linsolve(A,b)

