

## scalar fun times vector

*{T<:Union(Number,IFun)}(f::IFun,v::Vector{T})=typeof(f)[f.*v[k] for k=1:length(v)]
*{T<:Union(Number,IFun)}(v::Vector{T},f::IFun)=typeof(f)[v[k].*f for k=1:length(v)]
*(f::IFun,v::Vector{Any})=typeof(f)[f.*v[k] for k=1:length(v)]
*(v::Vector{Any},f::IFun)=typeof(f)[v[k].*f for k=1:length(v)]


## Vector of fun routines

function coefficients{N<:Number}(f::Vector{IFun{N}},m...)
    n=mapreduce(length,max,f)
    R=zeros(n,length(f))
    for k=1:length(f)
        R[1:length(f[k]),k]=coefficients(f,m)
    end
    R
end


function coefficients{T<:FFun}(B::Vector{T})
    m=mapreduce(length,max,B)
    fi=mapreduce(f->firstindex(f.coefficients),min,B)

    n=length(B)
    ret = zeros(Complex{Float64},m,length(B))
    for j=1:n
        for k=firstindex(B[j].coefficients):lastindex(B[j].coefficients)
            ret[k - fi + 1,j] = B[j].coefficients[k]
        end
    end
  
    ret
end


#TODO: Fun*vec should be Array[IFun]
*{T<:IFun}(v::Vector{T},a::Vector)=IFun(coefficients(v)*a,first(v).domain) 
function *{T<:FFun}(v::Vector{T},a::Vector)
    fi=mapreduce(f->firstindex(f.coefficients),min,v)
    FFun(ShiftVector(coefficients(v)*a,1-fi),first(v).domain) 
end

## evaluation


#TODO: fix for complex 
evaluate{T<:AbstractFun}(A::Vector{T},x::Real)=Float64[real(A[k][x]) for k=1:length(A)]


function evaluate{T<:IFun}(A::Vector{T},x::Vector{Float64})
    x = tocanonical(first(A),x)

    n=length(x)
    ret=Array(Float64,length(A),n)
    
    bk=Array(Float64,n)
    bk1=Array(Float64,n)
    bk2=Array(Float64,n)
    
    for k=1:length(A)
        bkr=clenshaw(A[k].coefficients,x,bk,bk1,bk2)
        
        for j=1:n
            ret[k,j]=bkr[j]
        end
    end
    
    ret
end

function evaluate{T<:FFun}(A::Vector{T},x::Vector{Float64})
    x = tocanonical(first(A),x)

    n=length(x)
    ret=Array(Float64,length(A),n)
    
    for k=1:length(A)
        bk=horner(A[k].coefficients,x)
        
        for j=1:n
            ret[k,j]=bk[j]
        end
    end
    
    ret
end

