


## Vector of fun routines

function coefficients{N,F}(::Type{N},f::Vector{F},o...)
    if isempty(f)
        return Array(N,0,0)
    end

    n=mapreduce(length,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[1:length(f[k]),k]=coefficients(f[k],o...)
    end
    R
end


scalarorfuntype{S,T<:Number}(::Fun{S,T})=T
scalarorfuntype{T<:Number}(::T)=T
scalarorfuntype{T<:Number}(b::Vector{T})=T
scalarorfuntype(b::Vector{Any})=promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::Vector{F})=promote_type(map(scalarorfuntype,b)...)


coefficients{F<:Fun}(Q::Vector{F},o...)=coefficients(scalarorfuntype(Q),Q,o...)
coefficients(Q::Vector{Any})=(@assert isempty(Q); zeros(0,0))

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


function values{D,N}(f::Vector{Fun{D,N}})
    n=mapreduce(length,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[:,k] = values(pad(f[k],n))
    end
    R
end

function values{D,T}(p::Array{Fun{D,T},2})
    @assert size(p)[1] == 1

   values(vec(p))
end







## evaluation


#TODO: fix for complex
evaluate{T<:Fun}(A::Vector{T},x::Number)=typeof(first(A)(x))[A[k](x) for k=1:length(A)]
evaluate{T<:Fun}(A::Array{T},x::Number)=typeof(first(A)(x))[A[k,j](x) for k=1:size(A,1),j=1:size(A,2)]


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



 for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
     @eval begin
         function ($op){T<:Number,V<:Number,D}(A::Array{T,2}, p::Vector{Fun{D,V}})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{D,promote_type(T,V)},size(cfs,1))
             for i = 1:size(cfs,1)
                 ret[i] = chop!(Fun(vec(cfs[i,:]),first(p).space),eps())
             end
             ret
         end

         function ($op){T<:Number,D}(p::Vector{Fun{D,T}},A::Array{T,2})
             cfs=$op(A,coefficients(p).')
             ret = Array(Fun{D,T},size(cfs,1))
             for i = 1:size(cfs,1)
                 ret[i] = chop!(Fun(vec(cfs[i,:]),first(p).space),eps())
             end
             ret
         end
     end
 end


#Allow vecfun + constvec, etc
#can't just promote constant vector to a vector-valued fun because don't know the domain.
for op = (:+,:-,:.*,:./)
    @eval begin
        ($op){T<:Number,S,V}(f::Fun{S,V},c::Array{T})=devec($op(vec(f),c))
        ($op){T<:Number,S,V}(c::Array{T},f::Fun{S,V})=devec($op(c,vec(f)))
    end
end
