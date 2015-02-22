## This file has methods to stably divide the singularity


# function dirichlet_transform{T<:Number}(s,v::Vector{T})
#     n=length(v)
#     w=zeros(T,n-1)
#     w[n-1]=v[n]
#     for k=n-2:-1:1
#         @inbounds w[k]=v[k+1] + s.*w[k+1] 
#     end
#     
#    @assert abs(s.*w[1] + v[1]) < 10000*n*eps()
#     
#     w
# end
# 
# function idirichlet_transform{T<:Number}(s,v::Vector{T})
#     n=length(v)
#     w=zeros(T,n+1)
#     w[1]=-s.*v[1]
#     for k=2:n
#         @inbounds w[k]=v[k-1] - s.*v[k] 
#     end
#     
#     w[n+1]=v[n]
#     
#     w
# end

# 
# function dirichletrange_divide_singularity{T<:Number}(s,v::Vector{T})
#     n=length(v)
#     w=zeros(T,n)
#     w[n]=-2s*v[n]
#     
#     for k=n-1:-1:2
#         @inbounds w[k]=-2s*(v[k] - .5w[k+1])
#     end
#     
#     w[1]=-s*v[1]
#     if n≥2
#         w[1]+=0.5s*w[2]
#     end
#     
#     w    
# end
# 
# function dirichlet_divide_singularity{T<:Number}(s,v::Vector{T})
#     n=length(v)
#     w=zeros(T,n-1)
#     w[n-1]=-2s*v[n]
#     w[n-2]=-2s*(v[n-1]-w[n-1])
#     
#     for k=n-3:-1:1
#         @inbounds w[k]=-2s*(v[k+1] - w[k+1] + .5s*w[k+2])
#     end
#     
#    @assert abs(1.5w[1] - .5w[2] - v[1]) < 10000n*eps()
#     
#     w    
# end


function dirichlet_divide_singularity{T<:Number}(b::Bool,v::Vector{T})
    n=length(v)
    w=zeros(T,n-1)
    s=b?1:-1
    w[n-1]=-2s*v[n]
    
    for k=n-2:-1:2
        @inbounds w[k]=-2s*(v[k+1] - .5w[k+1])
    end
    
    @inbounds w[1]=-s*(v[2] - .5w[2])    
    
    
    w    
end



#divide_singularity(s,v::Vector)=dirichletrange_divide_singularity(s,dirichlet_transform(s,v))
#divide_singularity(s,v::Vector)=idirichlet_transform(s,dirichlet_transform(s,dirichletrange_divide_singularity(s,dirichlet_transform(s,v))))
#divide_singularity(s,f::Fun)=Fun(divide_singularity(s,f.coefficients),f.space)




# 
# 
# ## both left and right dirichlet_transform
# converts to f_0 T_0 + f_1 T_1 +  \sum f_k (T_k - T_{k-2})
# function dirichlettransform!(w::Vector)
#     for k=length(w)-2:-1:1
#         @inbounds w[k] += w[k+2] 
#     end
#     
#     w
# end
# 
# function idirichlettransform!(w::Vector)
#     for k=3:length(w)
#         @inbounds w[k-2]-= w[k] 
#     end
#     
#     w
# end
# 
# dirichlettransform(v::Vector)=dirichlettransform!(copy(v))
# idirichlettransform(v::Vector)=idirichlettransform!(copy(v))



# dirichletrange_divide_singularity divides by (1-x^2)
# where the coefficients are ChebyshevDirichlet{1,1}
# assume that v[1] == v[2] ==0.
# the output coefficients are Chebyshev


function dirichlet_divide_singularity{T<:Number}(v::Vector{T})
    n=length(v)
    w=zeros(T,n-2)
    w[n-2]=-4v[n]
    if n>3
        w[n-3]=-4v[n-1]
    end
    
    for k=n-4:-1:2
        @inbounds w[k]=-4*(v[k+2] - .25w[k+2])
    end

    w[1]=-2v[3]
    if n≥5
        w[1]+=0.5w[3]
    end
    
    w    
end


divide_singularity(s::Bool,v::Vector)=dirichlet_divide_singularity(s,spaceconversion(v,Chebyshev,s?ChebyshevDirichlet{0,1}:ChebyshevDirichlet{1,0}))
divide_singularity(s::Integer,v::Vector)=divide_singularity(s==1,v)

divide_singularity(v::Vector)=dirichlet_divide_singularity(spaceconversion(v,Chebyshev,ChebyshevDirichlet{1,1}))


divide_singularity{T}(f::Fun{Chebyshev,T})=Fun(divide_singularity(f.coefficients),f.space)
divide_singularity{T}(s,f::Fun{Chebyshev,T})=Fun(divide_singularity(s,f.coefficients),f.space)
