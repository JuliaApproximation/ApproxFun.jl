

export PartialInverseOperator


struct PartialInverseOperator{T<:Number,CO<:CachedOperator,BI} <: Operator{T}
    cache::CO
    bandwidths::BI
end

function PartialInverseOperator(CO::CachedOperator{T},bandwidths) where T<:Number
    @assert istriu(CO) # || istril(CO)
    return PartialInverseOperator{T,typeof(CO),typeof(bandwidths)}(CO,bandwidths)
end

PartialInverseOperator(B::Operator, bandwidths) = PartialInverseOperator(cache(B), bandwidths)
PartialInverseOperator(B::Operator) = PartialInverseOperator(B, bandwidths(B))

convert(::Type{Operator{T}},A::PartialInverseOperator) where {T}=PartialInverseOperator(convert(Operator{T},A.cache), A.bandwidths)

domainspace(P::PartialInverseOperator)=rangespace(P.cache)
rangespace(P::PartialInverseOperator)=domainspace(P.cache)
domain(P::PartialInverseOperator)=domain(domainspace(P))
bandwidths(P::PartialInverseOperator) = P.bandwidths

function getindex(P::PartialInverseOperator,k::Integer,j::Integer)
    b = bandwidth(P.cache, 2)
    if k == j
        return inv(P.cache[k,k])
    elseif j > k
        t = zero(T)
        for i = max(k,j-b-1):j-1
            t += ret[k,i]*P.cache[i,j]
        end
        return -t/P.cache[j,j]
    else
        return zero(eltype(P))
    end
end
