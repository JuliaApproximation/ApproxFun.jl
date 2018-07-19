##Operators

for op in (:Derivative,:Integral)
    @eval begin
        function ($op)(d::AbstractVector{T}) where T<:IntervalDomain
            n=length(d)
            R=zeros(Operator{mapreduce(eltype,promote_type,d)},n,n)
            for k=1:n
                R[k,k]=$op(d[k])
            end

            R
        end
    end
end

function Evaluation(d::AbstractVector{T},x...) where T<:IntervalDomain
    n=length(d)
    R=zeros(Operator{mapreduce(eltype,promote_type,d)},n,n)
    for k=1:n
        R[k,k]=Evaluation(d[k],x...)
    end

    R
end


## Construction
function diagm_container(kv::Pair{<:Integer,<:AbstractVector{O}}...) where O<:Operator
    T = mapreduce(x -> mapreduce(eltype,promote_type,x.second),
                  promote_type, kv)
    n = mapreduce(x -> length(x.second) + abs(x.first), max, kv)
    zeros(Operator{T}, n, n)
end

##TODO: unify with other blockdiag
function blockdiag(d1::AbstractVector{T},d2::AbstractVector{T}) where T<:Operator
  if isempty(d1)&&isempty(d2)
    error("Empty blockdiag")
  end
  if isempty(d1)
    TT=mapreduce(eltype,promote_type,d2)
  elseif isempty(d2)
    TT=mapreduce(eltype,promote_type,d1)
  else
     TT=promote_type(mapreduce(eltype,promote_type,d1),
                    mapreduce(eltype,promote_type,d2))
  end

  D=zeros(Operator{TT},length(d1)+length(d2),2)
  D[1:length(d1),1]=d1
  D[length(d1)+1:end,2]=d2
  D
end

blockdiag(a::Operator,b::Operator) = blockdiag(Operator{promote_type(eltype(a),eltype(b))}[a],
                                                Operator{promote_type(eltype(a),eltype(b))}[b])

## broadcase

broadcast(::typeof(*),A::AbstractArray{N},D::Operator) where {N<:Number} =
    Operator{promote_type(N,eltype(D))}[A[k,j]*D for k=1:size(A,1),j=1:size(A,2)]
broadcast(::typeof(*),D::Operator,A::AbstractArray{N}) where {N<:Number}=A.*D
