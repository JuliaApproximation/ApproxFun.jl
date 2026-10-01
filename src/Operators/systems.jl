##Operators

for op in (:Derivative,:Integral)
    @eval begin
        function ($op){T<:IntervalDomain}(d::Vector{T})
            n=length(d)
            R=zeros(Operator{mapreduce(eltype,promote_type,d)},n,n)
            for k=1:n
                R[k,k]=$op(d[k])
            end

            R
        end
    end
end




function Evaluation{T<:IntervalDomain}(d::Vector{T},x...)
    n=length(d)
    R=zeros(Operator{mapreduce(eltype,promote_type,d)},n,n)
    for k=1:n
        R[k,k]=Evaluation(d[k],x...)
    end

    R
end



## Construction

function Base.diagm{T<:Operator}(d::Vector{T})
    D=zeros(Operator{mapreduce(eltype,promote_type,d)},length(d),length(d))
    for k=1:length(d)
        D[k,k]=d[k]
    end
    D
end

##TODO: unify with other blkdiag
function Base.blkdiag{T<:Operator}(d1::Vector{T},d2::Vector{T})
  if isempty(d1)&&isempty(d2)
    error("Empty blkdiag")
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

Base.blkdiag(a::Operator,b::Operator) = blkdiag(Operator{promote_type(eltype(a),eltype(b))}[a],
                                                Operator{promote_type(eltype(a),eltype(b))}[b])

## broadcase

.*{N<:Number}(A::Array{N},D::Operator)=Operator{promote_type(N,eltype(D))}[A[k,j]*D for k=1:size(A,1),j=1:size(A,2)]
.*{N<:Number}(D::Operator,A::Array{N})=A.*D
