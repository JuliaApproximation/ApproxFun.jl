
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
points(d::DomainSpace,n)=points(domain(d),n)

##Check domain compatibility


##Check domain compatibility

domainscompatible(a::DomainSpace,b::DomainSpace) = a.domain == AnyDomain() || b.domain == AnyDomain() || a.domain == b.domain


domain(A::DomainSpace)=A.domain # assume it has a field domain

##max space

maxspace(a::AnySpace,b::AnySpace)=a
maxspace(a::FunctionSpace,b::AnySpace)=a
maxspace(b::AnySpace,a::FunctionSpace)=a
function maxspace(a::FunctionSpace,b::FunctionSpace)
    @assert a==b    
    a
end

minspace(a::AnySpace,b::AnySpace)=a
minspace(a::FunctionSpace,b::AnySpace)=a
minspace(b::AnySpace,a::FunctionSpace)=a
function minspace(a::FunctionSpace,b::FunctionSpace)
    @assert a==b
    a
end


##TODO: Do we need both max and min?
function findmindomainspace(ops::Vector)
    sp = AnySpace()
    
    for op in ops
        sp = minspace(sp,domainspace(op))
    end
    
    sp
end

function findmaxrangespace(ops::Vector)
    sp = AnySpace()
    
    for op in ops
        sp = maxspace(sp,rangespace(op))
    end
    
    sp
end







