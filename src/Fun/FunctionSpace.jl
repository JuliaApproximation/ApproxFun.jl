
abstract FunctionSpace

##TODO: Confusing: two uses of domain

 # T tells whether the basis is real (cos/sin) or complex
abstract DomainSpace{T} <: FunctionSpace
abstract IntervalDomainSpace <: DomainSpace{Float64}     # We assume basis is real
abstract PeriodicDomainSpace{T} <: DomainSpace{T}       



export FunctionSpace, ChebyshevSpace, domainspace, rangespace, maxspace, minspace,Space


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




## TODO: remove zeros
Base.zero(S::FunctionSpace)=zeros(S)  
Base.zero{T<:Number}(::Type{T},S::FunctionSpace)=zeros(T,S)
Base.zeros{T<:Number}(::Type{T},S::FunctionSpace)=Fun(zeros(T,1),S)
Base.zeros(S::FunctionSpace)=Fun(zeros(1),S)




## Finite dimensional spaces



immutable VectorSpace{d} <: FunctionSpace
end

typealias ScalarSpace VectorSpace{1}

=={d}(::VectorSpace{d},::VectorSpace{d})=true
spacescompatible{d}(::VectorSpace{d},::VectorSpace{d})=true