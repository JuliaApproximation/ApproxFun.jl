


## Vector of fun routines

function coefficients{N,D}(f::Vector{Fun{N,D}},o...)
    n=mapreduce(length,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[1:length(f[k]),k]=coefficients(f[k],o...)
    end
    R
end




# function coefficients{T<:FFun}(B::Vector{T})
#     m=mapreduce(length,max,B)
#     fi=mapreduce(f->firstindex(f.coefficients),min,B)
# 
#     n=length(B)
#     ret = zeros(Complex{Float64},m,length(B))
#     for j=1:n
#         for k=firstindex(B[j].coefficients):lastindex(B[j].coefficients)
#             ret[k - fi + 1,j] = B[j].coefficients[k]
#         end
#     end
#   
#     ret
# end


function values{N}(f::Vector{Fun{N}})
    n=mapreduce(length,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[:,k] = values(pad(f[k],n))
    end
    R
end

function values{T}(p::Array{Fun{T},2})
    @assert size(p)[1] == 1

   values(vec(p))
end







## evaluation


#TODO: fix for complex 
evaluate{T<:Fun}(A::Vector{T},x::Real)=Float64[real(A[k][x]) for k=1:length(A)]
evaluate{T<:Fun}(A::Array{T},x::Real)=Float64[real(A[k,j][x]) for k=1:size(A,1),j=1:size(A,2)]


function evaluate{T<:Fun}(A::Vector{T},x::Vector{Float64})
    x = tocanonical(first(A),x)

    n=length(x)
    ret=Array(Float64,length(A),n)
    
    cplan=ClenshawPlan(Float64,n)
    
    for k=1:length(A)
        bkr=clenshaw(A[k].coefficients,x,cplan)
        
        for j=1:n
            ret[k,j]=bkr[j]
        end
    end
    
    ret
end

# function evaluate{T<:FFun}(A::Vector{T},x::Vector{Float64})
#     x = tocanonical(first(A),x)
# 
#     n=length(x)
#     ret=Array(Float64,length(A),n)
#     
#     for k=1:length(A)
#         bk=horner(A[k].coefficients,x)
#         
#         for j=1:n
#             ret[k,j]=bk[j]
#         end
#     end
#     
#     ret
# end



## Algebra

## scalar fun times vector

# *{T<:Union(Number,Fun)}(f::Fun,v::Vector{T})=typeof(f)[f.*v[k] for k=1:length(v)]
# *{T<:Union(Number,Fun)}(v::Vector{T},f::Fun)=typeof(f)[v[k].*f for k=1:length(v)]
# *(f::Fun,v::Vector{Any})=typeof(f)[f.*v[k] for k=1:length(v)]
# *(v::Vector{Any},f::Fun)=typeof(f)[v[k].*f for k=1:length(v)]
#  

#*{T<:Fun}(v::Vector{T},a::Vector)=Fun(coefficients(v)*a,first(v).space) 

# 
# function *{T<:FFun}(v::Vector{T},a::Vector)
#     fi=mapreduce(f->firstindex(f.coefficients),min,v)
#     FFun(ShiftVector(coefficients(v)*a,1-fi),first(v).domain) 
# end

# *{T<:Fun}(v::Vector{T},a::Vector)=Fun(coefficients(v)*a,first(v).space) 
# *{T<:Fun}(v::Vector{T},a::Number)=T[vk*a for vk in v]
# *{T<:Fun}(a::Number,v::Vector{T})=T[vk*a for vk in v]
# 
# ## Need to catch A*p, A'*p, A.'*p
# ##TODO: A may not be same type as p


dotu{N<:Real,D,T}(f::Vector{Fun{N,D}},g::Vector{T})=dot(f,g)
dotu{D,T}(f::Vector{Fun{Complex{Float64},D}},g::Vector{T})=dot(conj(f),g)

 for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
     @eval begin
         function ($op){T<:Number,D}(A::Array{T,2}, p::Vector{Fun{T,D}})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{T,D},size(cfs,1))
             for i = 1:size(A)[1]
                 ret[i] = Fun(vec(cfs[i,:]),p[i].space)
             end
             ret    
         end
         function ($op){T<:Number,V<:Number,D}(A::Array{T,2}, p::Vector{Fun{V,D}})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{promote_type(T,V),D},size(cfs,1))
             for i = 1:size(A)[1]
                 ret[i] = Fun(vec(cfs[i,:]),p[i].space)
             end
             ret    
         end      
         
         function ($op){T<:Number,D}(p::Vector{Fun{T,D}},A::Array{T,2})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{T,D},size(cfs,1))
             for i = 1:size(A)[1]
                 ret[i] = Fun(vec(cfs[i,:]),p[i].space)
             end
             ret    
         end
         function ($op){T<:Number,V<:Number,D}(A::Array{T,2}, p::Vector{Fun{V,D}})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{promote_type(T,V),D},size(cfs,1))
             for i = 1:size(A)[1]
                 ret[i] = Fun(vec(cfs[i,:]),p[i].space)
             end
             ret    
         end                
     end
 end