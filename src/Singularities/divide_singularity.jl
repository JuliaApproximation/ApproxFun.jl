## This file has methods to stably divide the singularity


function dirichlet_transform{T<:Number}(s,v::Vector{T})
    n=length(v)
    w=zeros(T,n-1)
    w[n-1]=v[n]
    for k=n-2:-1:1
        @inbounds w[k]=v[k+1] + s.*w[k+1] 
    end
    
#    @assert abs(s.*w[1] + v[1]) < 10000*n*eps()
    
    w
end

function idirichlet_transform{T<:Number}(s,v::Vector{T})
    n=length(v)
    w=zeros(T,n+1)
    w[1]=-s.*v[1]
    for k=2:n
        @inbounds w[k]=v[k-1] - s.*v[k] 
    end
    
    w[n+1]=v[n]
    
    w
end

function dirichlet_divide_singularity{T<:Number}(s,v::Vector{T})
    n=length(v)
    w=zeros(T,n-1)
    w[n-1]=-2s*v[n]
    w[n-2]=-2s*(v[n-1]-w[n-1])
    
    for k=n-3:-1:1
        @inbounds w[k]=-2s*(v[k+1] - w[k+1] + .5s*w[k+2])
    end
    
#    @assert abs(1.5w[1] - .5w[2] - v[1]) < 10000n*eps()
    
    w    
end

divide_singularity(s,v::Vector)=idirichlet_transform(s,dirichlet_divide_singularity(s,dirichlet_transform(s,v)))
divide_singularity(s,f::IFun)=IFun(divide_singularity(s,f.coefficients),f.domain)
divide_singularity(f::IFun)=divide_singularity(+1,divide_singularity(-1,f))




## both left and right dirichlet_transform

function dirichlet_transform{T<:Number}(v::Vector{T})
    n=length(v)
    w=zeros(T,n-2)
    w[n-2]=v[n]
    w[n-3]=v[n-1]    
    for k=n-4:-1:1
        @inbounds w[k]=v[k+2] + w[k+2] 
    end
    
    w
end



function idirichlet_transform{T<:Number}(v::Vector{T})
    n=length(v)
    w=zeros(T,n+2)
    w[1]=-v[1]
    w[2]=-v[2]    
    for k=3:n
        @inbounds w[k]=v[k-2] - v[k] 
    end
    
    w[n+1]=v[n-1]
    w[n+2]=v[n]
    
    w
end
