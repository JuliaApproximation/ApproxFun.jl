
abstract OperatorSpace

## Any is allowed as a Space
typealias Space Union(OperatorSpace,DataType)


type VectorSpace <: OperatorSpace
    dimension::Int
end

==(a::VectorSpace,b::VectorSpace)= a.dimension==b.dimension

##Check domain compatibility

domainscompatible(a::DataType,b::DataType) = a == Any && b==Any
domainscompatible(a::DataType,b) = a == Any
domainscompatible(a,b::DataType) = b == Any
domainscompatible(a,b) = a==b


spacescompatible(a::DataType,b::DataType) = a == Any && b==Any
spacescompatible(a::DataType,b) = a == Any
spacescompatible(a,b::DataType) = b == Any
spacescompatible(a,b) = a==b

##Default is Any

rangespace(A::Operator)=Any
domainspace(A::Operator)=Any

##max space


Base.max(a::OperatorSpace,b::DataType)=a
Base.max(b::DataType,a)=a
function Base.max(a::OperatorSpace,b::OperatorSpace)
    @assert a==b
    
    a
end


Base.min(a::OperatorSpace,b::DataType)=a
Base.min(b::DataType,a)=a
function Base.min(a::OperatorSpace,b::OperatorSpace)
    @assert a==b
    
    a
end

function findmindomainspace(ops::Vector)
    sp = Any
    
    for op in ops
        sp = min(sp,domainspace(op))
    end
    
    sp
end

function findmaxrangespace(ops::Vector)
    sp = Any
    
    for op in ops
        sp = max(sp,rangespace(op))
    end
    
    sp
end



function promoterangespace{T<:Operator}(ops::Vector{T})
    k=findmaxrangespace(ops)
    Operator[promoterangespace(op,k) for op in ops]
end




function promotedomainspace{T<:Operator}(ops::Vector{T})
    k=findmindomainspace(ops)
    T[promotedomainspace(op,k) for op in ops]
end

#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))

