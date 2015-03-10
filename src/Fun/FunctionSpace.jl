

export FunctionSpace, domainspace, rangespace, maxspace,Space,conversion_type


## 
# "enum" types used for whether the basis is
# Real or Complex.  AnyBasis is used when the answer
# is unknown.  
##

immutable RealBasis end
immutable ComplexBasis end
immutable AnyBasis end


promote_rule(::Type{RealBasis},::Type{ComplexBasis})=ComplexBasis
promote_rule(::Type{ComplexBasis},::Type{AnyBasis})=AnyBasis
promote_rule(::Type{RealBasis},::Type{AnyBasis})=AnyBasis

# coefficient_type(basis,valuetype) gives the type for coefficients
# for basis of type RealBasis/ComplexBasis and valuetype
# giving the type of function values
coefficient_type{T<:Complex}(::Type{ComplexBasis},::Type{T})=T
coefficient_type{T<:Real}(::Type{ComplexBasis},::Type{T})=Complex{T}
coefficient_type{T}(::Type{RealBasis},::Type{T})=T


#
# eltype for RealBasis/ComplexBasis gives the
# default type.  Maybe should be defaulteltype?
#

Base.eltype(::RealBasis)=Float64
Base.eltype(::ComplexBasis)=Complex{Float64}
Base.eltype(::AnyBasis)=Number

Base.eltype(::Type{RealBasis})=Float64
Base.eltype(::Type{ComplexBasis})=Complex{Float64}
Base.eltype(::Type{AnyBasis})=Number




# T is either RealBasis (cos/sin/polynomial) or ComplexBasis (laurent)
# D tells what canonical domain is (Interval/PeriodicInterval)
abstract FunctionSpace{T,D} #TODO should be able to write D<:Domain

typealias RealSpace{D} FunctionSpace{RealBasis,D}
typealias ComplexSpace{D} FunctionSpace{ComplexBasis,D}


Base.eltype{S}(::FunctionSpace{S})=eltype(S)
coefficient_type{S}(::FunctionSpace{S},T)=coefficient_type(S,T)

domain{T,D<:AnyDomain}(::FunctionSpace{T,D})=AnyDomain()




abstract AmbiguousSpace <: FunctionSpace{RealBasis,AnyDomain}


# AnySpace dictates that an operator can act on any space
# UnsetSpace dictates that an operator is not defined until 
#   its domainspace is promoted
# NoSpace is used to indicate no space exists for, e.g.,
# conversion_type

immutable AnySpace <: AmbiguousSpace end
immutable UnsetSpace <: AmbiguousSpace end
immutable NoSpace <: AmbiguousSpace end



#TODO: should it default to canonicalspace?
points(d::FunctionSpace,n)=points(domain(d),n)


##Check domain compatibility

Base.isapprox(a::Domain,b::Domain)=a==b
domainscompatible(a,b) = isa(domain(a),AnyDomain) || isa(domain(b),AnyDomain) || isapprox(domain(a),domain(b))

# Check whether spaces are the same, override when you need to check parameters
# This is used in place of == to support AnyDomain
spacescompatible{D<:FunctionSpace}(f::D,g::D)=error("Override spacescompatible for "*string(D))
spacescompatible(::AnySpace,::AnySpace)=true
spacescompatible(::UnsetSpace,::UnsetSpace)=true
spacescompatible(::NoSpace,::NoSpace)=true
spacescompatible(f,g)=false
==(A::FunctionSpace,B::FunctionSpace)=spacescompatible(A,B)&&domain(A)==domain(B)


# check a list of spaces for compatibility
function spacescompatible{T<:FunctionSpace}(v::Vector{T})
    for k=1:length(v)-1 
        if !spacescompatible(v[k],v[k+1])
            return false
        end
    end
    true
end

spacescompatible{T<:FunctionSpace}(v::Array{T})=spacescompatible(vec(v))



domain(A::FunctionSpace)=A.domain # assume it has a field domain



for op in (:tocanonical,:fromcanonical,:tocanonicalD,:fromcanonicalD)
    @eval ($op)(sp::FunctionSpace,x)=$op(domain(sp),x)
end




conversion_rule(a::FunctionSpace,b::FunctionSpace)=NoSpace()
function conversion_rule{S<:FunctionSpace}(a::S,b::S)
    if spacescompatible(a,b)
        a
    else
        NoSpace()
    end
end 

for FUNC in (:conversion_type,:maxspace)
    @eval begin
        $FUNC(::AnySpace,::UnsetSpace)=UnsetSpace()
        $FUNC(::UnsetSpace,::AnySpace)=UnsetSpace()        
    end

    for TYP in (:AnySpace,:UnsetSpace)
        @eval begin
            $FUNC(a::$TYP,b::$TYP)=a
            $FUNC(a::$TYP,b::FunctionSpace)=b
            $FUNC(a::FunctionSpace,b::$TYP)=a
        end
    end
    @eval $FUNC(b::AnySpace,a::FunctionSpace)=a
end


# gives a space c that has a banded conversion operator TO a and b
function conversion_type(a,b)
    if spacescompatible(a,b)
        a
    else
        cr=conversion_rule(a,b)
        cr==NoSpace()?conversion_rule(b,a):cr
    end
end







# gives a space c that has a banded conversion operator FROM a and b
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
#       coefficients(v::Vector,a,b) 
# converts from space a to space b
#       coefficients(v::Fun,a)
# is equivalent to coefficients(v.coefficients,v.space,a)
#       coefficients(v::Vector,a,b,c)
# uses an intermediate space b

coefficients(f::Vector,sp1::FunctionSpace,sp2::FunctionSpace,sp3::FunctionSpace)=coefficients(coefficients(f,sp1,sp2),sp2,sp3)

coefficients{T1<:FunctionSpace,T2<:FunctionSpace}(f::Vector,::Type{T1},::Type{T2})=coefficients(f,T1(),T2())
coefficients{T1<:FunctionSpace}(f::Vector,::Type{T1},sp2::FunctionSpace)=coefficients(f,T1(),sp2)
coefficients{T2<:FunctionSpace}(f::Vector,sp1::FunctionSpace,::Type{T2})=coefficients(f,sp1,T2())

## coefficients defaults to calling Conversion, otherwise it tries to pipe through Chebyshev


function coefficients{A<:FunctionSpace,B<:FunctionSpace}(f::Vector,a::A,b::B)
    ct=conversion_type(a,b) # gives a space that has a banded conversion to both a and b

    if spacescompatible(a,b)
        f
    elseif spacescompatible(ct,a)
        Conversion(a,b)*f  #TODO: Make * and \ consistent in return type
    elseif spacescompatible(ct,b)
        (Conversion(b,a)\f).coefficients
    else
        csp=canonicalspace(a)
        
        if spacescompatible(a,csp)# a is csp, so try b
            csp=canonicalspace(b)  
        end    
        if spacescompatible(a,csp)||spacescompatible(b,csp)# b is csp too, so we are stuck, try Fun constructor
            coefficients(Fun(x->Fun(f,a)[x],b))
        else
            coefficients(f,a,csp,b)
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







## rand

Base.rand(d::FunctionSpace,k...)=rand(domain(d),k...)
checkpoints(d::FunctionSpace)=checkpoints(domain(d))



## default transforms

function itransform{T}(S::FunctionSpace{T},cfs)
    csp=canonicalspace(S)
    if S==csp
        error("Override itransform(::"*string(typeof(S))*",cfs)")
    end    
    
    itransform(csp,coefficients(cfs,S,csp))
end

function transform{T}(S::FunctionSpace{T},vals)
    csp=canonicalspace(S)
    if S==csp
        error("Override transform(::"*string(typeof(S))*",vals)")
    end
    
    coefficients(transform(csp,vals),csp,S)
end


