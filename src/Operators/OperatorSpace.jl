
abstract OperatorSpace

export OperatorSpace, domainspace, rangespace

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



##Default is Any

rangespace(A::Operator)=Any
domainspace(A::Operator)=Any
domain(A::OperatorSpace)=A.domain # assume it has a field domain

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



## promotion


## Space Operator is used to wrap an Any space operator 
type SpaceOperator{T<:Number,O<:Operator{T},S<:OperatorSpace} <: BandedOperator{T}
    op::O
    space::S
#     
#     function SpaceOperator{T,O,S}(o::O,s::S)
#         @assert domainspace(o)==rangespace(o)==Any
#         new(o,s)
#     end
end

SpaceOperator{T<:Number,S<:OperatorSpace}(o::Operator{T},s::S)=SpaceOperator{T,typeof(o),S}(o,s)

domainspace(S::SpaceOperator)=S.space
rangespace(S::SpaceOperator)=S.space
addentries!(S::SpaceOperator,A,kr)=addentries!(S.op,A,kr)
bandinds(S::SpaceOperator)=bandinds(S.op)
domain(S::SpaceOperator)=domain(S.space)


for op in (:promoterangespace,:promotedomainspace)
    @eval begin
        ($op)(P::Operator,::Type{Any})=P
        ($op)(P::Operator,sp::OperatorSpace,::Type{Any})=SpaceOperator(P,sp)
    end
end

promoterangespace(P::Operator,sp::OperatorSpace)=promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::OperatorSpace)=promotedomainspace(P,sp,domainspace(P))
        
        
promoterangespace(P::Operator,sp::OperatorSpace,cursp::OperatorSpace)=(sp==cursp)?P:TimesOperator(ConversionOperator(cursp,sp),P)
promotedomainspace(P::Operator,sp::OperatorSpace,cursp::OperatorSpace)=(sp==cursp)?P:TimesOperator(P,ConversionOperator(sp,cursp))




function promoterangespace{T<:Operator}(ops::Vector{T})
    k=findmaxrangespace(ops)
    Operator[promoterangespace(op,k) for op in ops]
end




function promotedomainspace{T<:Operator}(ops::Vector{T})
    k=findmindomainspace(ops)
    Operator[promotedomainspace(op,k) for op in ops]
end

#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))





