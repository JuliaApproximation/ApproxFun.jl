
abstract FunctionSpace{T}

##TODO: Confusing: two uses of domain

 # T tells whether the basis is real (cos/sin) or complex
abstract DomainSpace{T} <: FunctionSpace{T}
abstract IntervalDomainSpace <: DomainSpace{Float64}     # We assume basis is real
abstract PeriodicDomainSpace{T} <: DomainSpace{T}       



export FunctionSpace, ChebyshevSpace, domainspace, rangespace, maxspace, minspace,Space


immutable ConstantSpace <: FunctionSpace{Float64}
end

domain(::ConstantSpace)=AnyDomain()


immutable AnySpace <: FunctionSpace{Float64}
end

immutable NoSpace <: FunctionSpace{Float64}
end

domain(::AnySpace)=AnyDomain()

#TODO: should it default to canonicalspace?
points(d::DomainSpace,n)=points(domain(d),n)



##Check domain compatibility

domainscompatible(a,b) = domain(a) == AnyDomain() || domain(b) == AnyDomain() || domain(a) == domain(b)

#Check whether spaces are the same, override when you need to check parameters
spacescompatible{D<:DomainSpace}(f::D,g::D)=domainscompatible(f,g) 
spacescompatible(f,g)=false
==(A::DomainSpace,B::DomainSpace)=spacescompatible(A,B)&&domain(A)==domain(B)


# check a list of spaces for compatibility
function spacescompatible{T<:FunctionSpace}(v::Vector{T})
    for k=1:length(v)-1 
        if !spacescompatible(v[k],v[k+1])
            return false
        end
    end
    true
end



domain(A::DomainSpace)=A.domain # assume it has a field domain

canonicaldomain{T<:IntervalDomainSpace}(::Type{T})=Interval()
canonicaldomain{T<:PeriodicDomainSpace}(::Type{T})=PeriodicInterval()




for op in (:tocanonical,:fromcanonical,:tocanonicalD,:fromcanonicalD)
    @eval ($op)(sp::DomainSpace,x)=$op(domain(sp),x)
end




conversion_rule(a::FunctionSpace,b::FunctionSpace)=NoSpace()
function conversion_rule{S<:FunctionSpace}(a::S,b::S)
    if spacescompatible(a,b)
        a
    else
        NoSpace()
    end
end 

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
        return a
    end
    
    cr=conversion_type(a,b)
    if cr==a
        return b
    elseif cr ==b
        return a
    end
    
    # check if its banded through canonicalspace
    cspa=canonicalspace(a)
    if cspa != a && maxspace(cspa,a)==cspa
        return maxspace(b,cspa)
    end
    
    cspb=canonicalspace(b)
    if cspb !=b && maxspace(cspb,b)==cspb
        return maxspace(a,cspb)            
    end
    
    NoSpace()
end


union_rule(a::FunctionSpace,b::FunctionSpace)=NoSpace()
function Base.union(a::FunctionSpace,b::FunctionSpace)
    if spacescompatible(a,b)
        return a
    end
    
    cr=union_rule(a,b)
    if cr==NoSpace()
        cr=union_rule(b,a)
    end
    
    
    if cr==NoSpace()
        cspa=canonicalspace(a)
        cspb=canonicalspace(b)
        if cspa!=a || cspb!=b
            cr=union(cspa,cspb)  #Max or min space?
        end
    end
    
    if cr==NoSpace()
        cr=maxspace(a,b)  #Max or min space?
    end    
    
    return cr
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
        
        if spacescompatible(a,csp)# a is csp, so try b
            csp=canonicalspace(b)  
        end    
        if spacescompatible(b,csp)# b is csp too, so we are stuck, try Fun constructor
            coefficients(Fun(x->Fun(f,a)[x],b))
        else
            spaceconversion(f,a,csp,b)
        end
    end
end




## TODO: remove zeros
Base.zero(S::FunctionSpace)=zeros(S)  
Base.zero{T<:Number}(::Type{T},S::FunctionSpace)=zeros(T,S)
Base.zeros{T<:Number}(::Type{T},S::FunctionSpace)=Fun(zeros(T,1),S)
Base.zeros(S::FunctionSpace)=Fun(zeros(1),S)

# catch all
Base.ones(S::FunctionSpace)=Fun(x->1.0,S)
Base.ones{T<:Number}(::Type{T},S::FunctionSpace)=Fun(x->one(T),S)
identity_fun(S::FunctionSpace)=Fun(x->x,S)


## Finite dimensional spaces



immutable VectorSpace{d} <: FunctionSpace{Float64}
end

typealias ScalarSpace VectorSpace{1}

=={d}(::VectorSpace{d},::VectorSpace{d})=true
spacescompatible{d}(::VectorSpace{d},::VectorSpace{d})=true




## rand

Base.rand(d::DomainSpace)=rand(domain(d))




## default transforms

function itransform{T}(S::FunctionSpace{T},cfs)
    csp=canonicalspace(S)
    itransform(csp,spaceconversion(cfs,S,csp))
end

function transform{T}(S::FunctionSpace{T},vals)
    csp=canonicalspace(S)
    spaceconversion(transform(csp,vals),csp,S)
end


