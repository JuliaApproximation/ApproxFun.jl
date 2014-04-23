
abstract OperatorSpace

## Any is allowed as a Space
typealias Space Union(OperatorSpace,DataType)


#Any is allowed as a domain
type UltrasphericalSpace{T<:Union(IntervalDomain,DataType)} <: OperatorSpace
    order::Int
    domain::T     
end

UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,Any)

==(a::UltrasphericalSpace,b::UltrasphericalSpace)=a.order==b.order && a.domain==b.domain

type FourierSpace{T<:Union(PeriodicDomain,DataType)} <: OperatorSpace
    domain::T
end

FourierSpace()=FourierSpace(Any)

==(a::FourierSpace,b::FourierSpace)= a.domain==b.domain

type VectorSpace <: OperatorSpace
    dimension::Int
end

==(a::VectorSpace,b::VectorSpace)= a.dimension==b.dimension

##Check domain compatibility

domainscompatible(a::UltrasphericalSpace,b::UltrasphericalSpace) = a.domain == Any || b.domain == Any || a.domain == b.domain
domainscompatible(a::DataType,b::DataType) = a == Any && b==Any
domainscompatible(a::DataType,b) = a == Any
domainscompatible(a,b::DataType) = b == Any
domainscompatible(a,b) = a==b

spacescompatible(a::UltrasphericalSpace,b::UltrasphericalSpace) = domainscompatible(a,b) && a.order >= b.order
spacescompatible(a::DataType,b::DataType) = a == Any && b==Any
spacescompatible(a::DataType,b) = a == Any
spacescompatible(a,b::DataType) = b == Any
spacescompatible(a,b) = a==b

##Default is Any

rangespace(A::Operator)=Any
domainspace(A::Operator)=Any

##max space



function Base.max(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order > b.order?a:b
end

Base.max(a::OperatorSpace,b::DataType)=a
Base.max(b::DataType,a)=a
function Base.max(a::OperatorSpace,b::OperatorSpace)
    @assert a==b
    
    a
end

function Base.min(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order < b.order?a:b
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


## Operator space manipulation

function promoterangespace(P::Operator,sp::Space)
    psp = rangespace(P)
    
    if psp == Any || sp == Any || psp == sp
        P
    else
        @assert typeof(psp) <: UltrasphericalSpace
        @assert typeof(sp) <: UltrasphericalSpace
        @assert spacescompatible(sp,psp)
        
        if psp.order == sp.order
            P
        else
            ConversionOperator(psp.order:sp.order)*P
        end
    end
end


function promoterangespace{T<:Operator}(ops::Vector{T})
    k=findmaxrangespace(ops)
    T[promoterangespace(op,k) for op in ops]
end



function promotedomainspace(P::Operator,sp::Space)
    psp = domainspace(P)
    
    if psp == Any || sp == Any || psp == sp
        P
    else
        @assert typeof(psp) <: UltrasphericalSpace
        @assert typeof(sp) <: UltrasphericalSpace
        @assert spacescompatible(psp,sp)
        
        
        if psp.order == sp.order
            P
        else
            P*ConversionOperator(sp.order:psp.order)
        end
    end
end



function promotedomainspace{T<:Operator}(ops::Vector{T})
    k=findmindomainspace(ops)
    T[promotedomainspace(op,k) for op in ops]
end

#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))

##TODO: remove?
promotespaces(op::BandedOperator,od::Range1)=promoterangespace(promotedomainspace(op,UltrasphericalSpace(od[1])),UltrasphericalSpace(od[end]))
