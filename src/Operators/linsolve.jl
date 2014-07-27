## Linear Solve


IFun_coefficients(b::Vector,sp)=vcat(map(f-> typeof(f)<: IFun? coefficients(f,sp) :  f,b)...)
FFun_coefficients(b::Vector)=vcat(map(f-> typeof(f)<: FFun? interlace(f.coefficients) :  interlace(f),b)...) #Assume only FFun or ShiftVector

function IFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    u=adaptiveqr(A,IFun_coefficients(b,rangespace(A[end]).order),tolerance,maxlength)  ##TODO: depends on ordering of A
    
    IFun(u,commondomain(A,b))
end

function FFun_linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    @assert length(A) == 1

    u=adaptiveqr([interlace(A[1])],FFun_coefficients(b),tolerance,maxlength)
    
    FFun(deinterlace(u),commondomain(A,b))    
end

function linsolve{T<:Operator}(A::Vector{T},b::Vector;tolerance=0.01eps(),maxlength=1000000)
    d=commondomain(A,b)

    if typeof(d) <: IntervalDomain
        IFun_linsolve(A,b;tolerance=tolerance,maxlength=maxlength)
    elseif typeof(d) <: PeriodicDomain
        FFun_linsolve(A,b;tolerance=tolerance,maxlength=maxlength)
    else
        adaptiveqr(A,b,tolerance,maxlength)
    end    
end


##Todo nxn operator
 function linsolve{T<:Operator,M<:Number}(A::Array{T,2},b::Vector{M};tolerance=0.01eps(),maxlength=1000000)
    m = size(A,2)
 
     ret=adaptiveqr(interlace(A),b,tolerance,maxlength)  #Given just an array, we don't know how to interlace
                                                         #so assume user knows
                                     
                                     
     IFun[IFun(ret[k:m:end],commondomain(A[:,k])) for k=1:m]
 end
 
function linsolve{T<:Operator}(A::Array{T,2},b::Vector{Any};kwds...)
    m,n=size(A)

    br=m-n

    l=mapreduce(length,max,b[br+1:end])

    r=zeros(Float64,br+n*l)  ##TODO: support complex
    
    r[1:br]=b[1:br]
    
    for k=br+1:m
        sp=findmaxrangespace([A[k,:]...]).order
        if typeof(b[k])<:AbstractFun
            r[k:n:end]=pad(coefficients(b[k],sp),l)
        else  #type is vector
            r[k:n:end]=pad(b[k],l)
        end
    end

    linsolve(A,r;kwds...)
end 


linsolve(A::Operator,b::Vector;kwds...)=linsolve([A],b;kwds...)
linsolve(A,b;kwds...)=linsolve(A,[b];kwds...)


\{T<:Operator}(A::Array{T,2},b::Vector)=linsolve(A,b)
\{T<:Operator}(A::Vector{T},b::Vector)=linsolve(A,b)
\(A::Operator,b)=linsolve(A,b)

