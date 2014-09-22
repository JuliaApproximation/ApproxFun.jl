
abstract FunctionSpace

##TODO: Confusing: two uses of domain
abstract DomainSpace <: FunctionSpace
abstract IntervalDomainSpace <: DomainSpace     # We assume basis is real
abstract PeriodicDomainSpace{T} <: DomainSpace  # T tells whether the basis is real (cos/sin) or complex


typealias RealDomainSpace Union(IntervalDomainSpace,PeriodicDomainSpace{Float64})


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

domainscompatible(a,b) = domain(a) == AnyDomain() || domain(b) == AnyDomain() || domain(a) == domain(b)

#Check whether spaces are the same, override when you need to check parameters
spacescompatible{D<:DomainSpace}(f::D,g::D)=domainscompatible(f,g) 
spacescompatible(f,g)=false
==(A::DomainSpace,B::DomainSpace)=spacescompatible(A,B)&&domain(A)==domain(B)

domain(A::DomainSpace)=A.domain # assume it has a field domain

canonicaldomain{T<:IntervalDomainSpace}(::Type{T})=Interval()
canonicaldomain{T<:PeriodicDomainSpace}(::Type{T})=PeriodicInterval()




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






## Conversion routines


## Space conversion default is through canonicalspace

spaceconversion(f::Vector,sp::FunctionSpace)=spaceconversion(f,canonicalspace(sp),sp)
spaceconversion(f::Vector,sp1::FunctionSpace,sp2::FunctionSpace,sp3::FunctionSpace)=spaceconversion(spaceconversion(f,sp1,sp2),sp2,sp3)


## spaceconversion defaults to calling Conversion, otherwise it tries to pipe through ChebyshevSpace

# function spaceconversion{A<:FunctionSpace}(f::Vector,a::A,b::A)
#     if spacescompatible(a,b)
#         f
# end

function spaceconversion{A<:FunctionSpace,B<:FunctionSpace}(f::Vector,a::A,b::B)
    ct=conversion_type(a,b)

    if spacescompatible(a,b)
        f
    elseif spacescompatible(ct,a)
        Conversion(a,b)*f  ##TODO: Make * and \ consistent in return type
    elseif spacescompatible(ct,b)
        (Conversion(b,a)\f).coefficients
    else
        csp=canonicalspace(a)
        if spacescompatible(a,csp)
            error("Override spaceconversion or implement Conversion from " * string(typeof(csp)) * " to " * string(B))
        elseif spacescompatible(b,csp)
            error("Override spaceconversion or implement Conversion from " * string(A) * " to " * string(typeof(csp)))
        else
            spaceconversion(f,a,csp,b)
        end
    end
end


