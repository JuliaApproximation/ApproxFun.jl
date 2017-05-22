


## Vector of fun routines

function coefficientmatrix{N,F}(::Type{N},f::Vector{F},o...)
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
scalarorfuntype{T<:Number}(b::Vector{T}) = T
scalarorfuntype(b::Vector{Any}) = promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::Vector{F}) = promote_type(map(scalarorfuntype,b)...)


coefficientmatrix{F<:Fun}(Q::Vector{F},o...)=coefficientmatrix(scalarorfuntype(Q),Q,o...)
coefficientmatrix(Q::Vector{Any})=(@assert isempty(Q); zeros(0,0))


function values(f::Vector{Fun{D,N,VN}}) where {D,N,VN}
    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[:,k] = values(pad(f[k],n))
    end
    R
end

function values(p::Array{Fun{D,T,VT},2}) where {D,T,VT}
    @assert size(p)[1] == 1

   values(vec(p))
end







## evaluation


#TODO: fix for complex
evaluate{T<:Fun}(A::AbstractArray{T},x::Number)=typeof(first(A)(x))[Akj(x) for Akj in A]



## Algebra

## scalar fun times vector


 for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
     @eval begin
         function ($op)(A::Array{T,2}, p::Vector{Fun{D,V,VT}}) where {T<:Number,V<:Number,D,VT}
             cfs=$op(A,coefficientmatrix(p).')
             ret = Vector{VFun{D,promote_type(T,V)}}(size(cfs,1))
             for i = 1:size(cfs,1)
                 ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
             end
             ret
         end

         function ($op)(p::Vector{Fun{D,T,VT}},A::Array{T,2}) where {T<:Number,D,VT}
             cfs=$op(A,coefficientmatrix(p).')
             ret = Vector{VFun{D,T}}(size(cfs,1))
             for i = 1:size(cfs,1)
                 ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
             end
             ret
         end
     end
 end


#Allow vecfun + constvec, etc
#can't just promote constant vector to a vector-valued fun because don't know the domain.
#TODO: Use Base?
for op = (:+,:-)
    @eval begin
        ($op)(f::Fun,c::AbstractArray{T}) where {T<:Number} = Fun($op(Array(f),c))
        ($op)(c::AbstractArray{T},f::Fun) where {T<:Number} = Fun($op(c,Array(f)))
    end
end
