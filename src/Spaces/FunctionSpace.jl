
abstract FunctionSpace

##TODO: Confusing: two uses of domain
abstract DomainSpace <: FunctionSpace
abstract IntervalDomainSpace <: DomainSpace
abstract PeriodicDomainSpace <: DomainSpace

export FunctionSpace, domainspace, rangespace, maxspace, minspace


type ConstantSpace <: FunctionSpace
end

domain(::ConstantSpace)=AnyDomain()

type AnySpace <: FunctionSpace
end

domain(::AnySpace)=AnyDomain()

## Any is allowed as a Space
typealias Space Union(FunctionSpace,DataType)



points(d::DomainSpace,n)=points(domain(d),n)

##Check domain compatibility

domainscompatible(a::DataType,b::DataType) = a == Any && b==Any
domainscompatible(a::DataType,b) = a == Any
domainscompatible(a,b::DataType) = b == Any
domainscompatible(a,b) = a==b  ##TODO: what is this suppose to do?

##Check domain compatibility

domainscompatible(a::DomainSpace,b::DomainSpace) = a.domain == Any || b.domain == Any || a.domain == b.domain


domain(A::DomainSpace)=A.domain # assume it has a field domain

##max space


maxspace(a::FunctionSpace,b::DataType)=a
maxspace(b::DataType,a)=a
function maxspace(a::FunctionSpace,b::FunctionSpace)
    @assert a==b
    
    a
end


minspace(a::FunctionSpace,b::DataType)=a
minspace(b::DataType,a)=a
function minspace(a::FunctionSpace,b::FunctionSpace)
    @assert a==b
    
    a
end


##TODO: Do we need both max and min?
function findmindomainspace(ops::Vector)
    sp = Any
    
    for op in ops
        sp = minspace(sp,domainspace(op))
    end
    
    sp
end

function findmaxrangespace(ops::Vector)
    sp = Any
    
    for op in ops
        sp = maxspace(sp,rangespace(op))
    end
    
    sp
end







