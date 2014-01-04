export coefficientmatrix

## scalar fun times vector

*{T<:Union(Number,IFun,Any)}(f::IFun,v::Vector{T})=typeof(f)[f.*v[k] for k=1:length(v)]
*{T<:Union(Number,IFun,Any)}(v::Vector{T},f::IFun)=typeof(f)[v[k].*f for k=1:length(v)]

## Vector of fun routines

function coefficientmatrix{T<:IFun}(B::Vector{T})
    m=mapreduce(length,max,B)
    n=length(B)
    ret = zeros(m,length(B))
    for j=1:n
        for k=1:length(B[j])
            ret[k,j] =  B[j].coefficients[k]
        end
    end
    
    ret
end

function coefficientmatrix{T<:FFun}(B::Vector{T})
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
*{T<:IFun}(v::Vector{T},a::Vector)=IFun(coefficientmatrix(v)*a,first(v).domain) 
function *{T<:FFun}(v::Vector{T},a::Vector)
    fi=mapreduce(f->firstindex(f.coefficients),min,v)
    FFun(ShiftVector(coefficientmatrix(v)*a,1-fi),first(v).domain) 
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
        bk=clenshaw(A[k].coefficients,x,bk,bk1,bk2)
        
        for j=1:n
            ret[k,j]=bk[j]
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

