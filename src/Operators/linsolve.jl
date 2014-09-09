## Linear Solve


IFun_coefficients(b::Vector,sp)=vcat(map(f-> isa(f,IFun)? coefficients(f,sp) :  f,b)...)
FFun_coefficients(b::Vector)=vcat(map(f-> isa(f,FFun)? interlace(f.coefficients) :  interlace(f),b)...) #Assume only FFun or ShiftVector

function IFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    u=adaptiveqr(A,IFun_coefficients(b,rangespace(A[end])),tolerance,maxlength)  ##TODO: depends on ordering of A
    
    IFun(u,commondomain(A,b))
end

function IFun_linsolve{T<:Operator,N<:Number}(A::Vector{T},b::Array{N,2};tolerance=0.01eps(),maxlength=1000000)
    u=adaptiveqr(A,b,tolerance,maxlength)  ##TODO: depends on ordering of A
    d=commondomain(A)
    IFun[IFun(u[:,k],d) for k=1:size(u,2)]
end



function FFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    @assert length(A) == 1

    u=adaptiveqr([interlace(A[1])],FFun_coefficients(b),tolerance,maxlength)
    
    FFun(deinterlace(u),commondomain(A,b))    
end

function linsolve{T<:Operator}(A::Vector{T},b::Array;tolerance=0.01eps(),maxlength=1000000)
    d=commondomain(A,b)

    if typeof(d) <: IntervalDomain
        IFun_linsolve(A,b;tolerance=tolerance,maxlength=maxlength)
    elseif typeof(d) <: PeriodicDomain
        FFun_linsolve(A,b;tolerance=tolerance,maxlength=maxlength)
    else
        adaptiveqr(A,b,tolerance,maxlength)
    end    
end


 function linsolve{T<:Operator,M<:Number}(A::Array{T,2},b::Vector{M};tolerance=0.01eps(),maxlength=1000000)
    m = size(A,2)
 
     ret=adaptiveqr(interlace(A),b,tolerance,maxlength)  #Given just an array, we don't know how to interlace
                                                         #so assume user knows, this is correct for bc rows
                                     
                                     
     IFun[IFun(ret[k:m:end],commondomain(A[:,k])) for k=1:m]
 end
 
 
function linsolve{T<:Operator,M<:Number}(A::Array{T,2},b::Array{M,2};tolerance=0.01eps(),maxlength=1000000)
    m = size(A,2)
    
    ret=adaptiveqr(interlace(A),b,tolerance,maxlength)  #Given just an array, we don't know how to interlace
                                                     #so assume user knows, this is correct for bc rows
                                 
                                 
    IFun[IFun(ret[k:m:end,j],commondomain(A[:,k])) for k=1:m,j=1:size(b,2)]
end
 

scalarorfuntype{T<:Number}(::IFun{T})=T
scalarorfuntype{T<:Number}(::T)=T
scalarorfuntype{T<:Number}(b::Vector{T})=T
scalarorfuntype(b::Vector{Any})=promote_type(map(scalarorfuntype,b)...)
 
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
        elseif isa(b[k],AbstractFun)
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

