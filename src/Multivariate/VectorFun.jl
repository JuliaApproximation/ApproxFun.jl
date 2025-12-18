


## Vector of fun routines

function coefficientmatrix{N,F}(::Type{N},f::Vector{F},o...)
    if isempty(f)
        return Array(N,0,0)
    end

    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[1:ncoefficients(f[k]),k]=coefficients(f[k],o...)
    end
    R
end


scalarorfuntype{S,T<:Number}(::Fun{S,T})=T
scalarorfuntype{T<:Number}(::T)=T
scalarorfuntype{T<:Number}(b::Vector{T})=T
scalarorfuntype(b::Vector{Any})=promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::Vector{F})=promote_type(map(scalarorfuntype,b)...)


coefficientmatrix{F<:Fun}(Q::Vector{F},o...)=coefficientmatrix(scalarorfuntype(Q),Q,o...)
coefficientmatrix(Q::Vector{Any})=(@assert isempty(Q); zeros(0,0))


function values{D,N}(f::Vector{Fun{D,N}})
    n=mapreduce(ncoefficients,max,f)
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
evaluate{T<:Fun}(A::AbstractArray{T},x::Number)=typeof(first(A)(x))[Akj(x) for Akj in A]


function evaluate{T<:Fun}(A::AbstractVector{T},x::AbstractVector)
    n=length(x)
    ret=Array(promote_type(eltype(x),mapreduce(eltype,promote_type,A)),length(A),n)

    for k=1:length(A)
        bkr=evaluate(A[k],x)

        for j=1:n
            ret[k,j]=bkr[j]
        end
    end

    ret
end



## Algebra

## scalar fun times vector


 for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
     @eval begin
         function ($op){T<:Number,V<:Number,D}(A::Array{T,2}, p::Vector{Fun{D,V}})
             cfs=$op(A,coefficientmatrix(p).')
             ret = Array(Fun{D,promote_type(T,V)},size(cfs,1))
             for i = 1:size(cfs,1)
                 ret[i] = chop!(Fun(vec(cfs[i,:]),first(p).space),eps())
             end
             ret
         end

         function ($op){T<:Number,D}(p::Vector{Fun{D,T}},A::Array{T,2})
             cfs=$op(A,coefficientmatrix(p).')
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
        ($op){T<:Number,S,V}(f::Fun{S,V},c::AbstractArray{T})=devec($op(vec(f),c))
        ($op){T<:Number,S,V}(c::AbstractArray{T},f::Fun{S,V})=devec($op(c,vec(f)))
    end
end
