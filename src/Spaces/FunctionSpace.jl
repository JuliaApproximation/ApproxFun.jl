
abstract FunctionSpace

##TODO: Confusing: two uses of domain
abstract DomainSpace <: FunctionSpace
abstract IntervalDomainSpace <: DomainSpace
abstract PeriodicDomainSpace <: DomainSpace

export FunctionSpace, ChebyshevSpace, domainspace, rangespace, maxspace, minspace


immutable ConstantSpace <: FunctionSpace
end

domain(::ConstantSpace)=AnyDomain()


immutable AnySpace <: FunctionSpace
end

immutable NoSpace <: FunctionSpace
end

domain(::AnySpace)=AnyDomain()
points(d::DomainSpace,n)=points(domain(d),n)



##Check domain compatibility

domainscompatible(a::DomainSpace,b::DomainSpace) = a.domain == AnyDomain() || b.domain == AnyDomain() || a.domain == b.domain


domain(A::DomainSpace)=A.domain # assume it has a field domain


for op in (:tocanonical,:fromcanonical,:tocanonicalD,:fromcanonicalD)
    @eval ($op)(sp::DomainSpace,x)=$op(domain(sp),x)
end




conversion_rule(a::FunctionSpace,b::FunctionSpace)=NoSpace()
conversion_rule{S<:FunctionSpace}(a::S,b::S)=a

function conversion_type(a,b)
    cr=conversion_rule(a,b)
    cr==NoSpace()?conversion_rule(b,a):cr
end


# gives a space c that has a banded conversion operator to a and b
minspace(a::AnySpace,b::AnySpace)=a
minspace(a::FunctionSpace,b::AnySpace)=a
minspace(b::AnySpace,a::FunctionSpace)=a
function minspace(a::FunctionSpace,b::FunctionSpace)
    if a==b
        a
    else
        conversion_type(a,b)
    end
end




# gives a space c that has a banded conversion operator from a and b
maxspace(a::AnySpace,b::AnySpace)=a
maxspace(a::FunctionSpace,b::AnySpace)=a
maxspace(b::AnySpace,a::FunctionSpace)=a
function maxspace(a::FunctionSpace,b::FunctionSpace)
    if a==b    
        a
    else
        cr=conversion_type(a,b)
        if cr==a
            b
        elseif cr ==b
            a
        else
            NoSpace()
        end
    end
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



