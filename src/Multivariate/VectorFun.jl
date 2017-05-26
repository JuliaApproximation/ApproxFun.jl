


## Vector of fun routines

function coefficientmatrix{N,F}(::Type{N},f::AbstractVector{F},o...)
    if isempty(f)
        return Matrix{N}(0,0)
    end

    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[1:ncoefficients(f[k]),k]=coefficients(f[k],o...)
    end
    R
end


scalarorfuntype{S,T<:Number}(::Fun{S,T}) = T
scalarorfuntype{T<:Number}(::T) = T
scalarorfuntype{T<:Number}(b::AbstractVector{T}) = T
scalarorfuntype(b::AbstractVector{Any}) = promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::AbstractVector{F}) = promote_type(map(scalarorfuntype,b)...)


coefficientmatrix{F<:Fun}(Q::AbstractVector{F},o...)=coefficientmatrix(scalarorfuntype(Q),Q,o...)
coefficientmatrix(Q::AbstractVector{Any})=(@assert isempty(Q); zeros(0,0))


function values(f::AbstractVector{Fun{D,N,VN}}) where {D,N,VN}
    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[:,k] = values(pad(f[k],n))
    end
    R
end

function values(p::AbstractMatrix{Fun{D,T,VT}}) where {D,T,VT}
    @assert size(p)[1] == 1

   values(vec(p))
end







## evaluation


#TODO: fix for complex
evaluate{T<:Fun}(A::AbstractArray{T},x::Number) =
    typeof(first(A)(x))[Akj(x) for Akj in A]



## Algebra
#
# ## scalar fun times vector
# #TODO: This is probably more efficient than the current implementation, but it
# # is too hacky to keep for now
#
#  for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
#      @eval begin
#          function ($op)(A::Array{T,2}, p::AbstractVector{Fun{D,V,VT}}) where {T<:Number,V<:Number,D,VT}
#              cfs=$op(A,coefficientmatrix(p).')
#              ret = Vector{VFun{D,promote_type(T,V)}}(size(cfs,1))
#              for i = 1:size(cfs,1)
#                  ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
#              end
#              ret
#          end
#
#          function ($op)(p::AbstractVector{Fun{D,T,VT}},A::Array{T,2}) where {T<:Number,D,VT}
#              cfs=$op(A,coefficientmatrix(p).')
#              ret = Vector{VFun{D,T}}(size(cfs,1))
#              for i = 1:size(cfs,1)
#                  ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
#              end
#              ret
#          end
#      end
#  end
