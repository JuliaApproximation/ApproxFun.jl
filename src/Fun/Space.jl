

export Space, domainspace, rangespace, maxspace,Space,conversion_type


##
# "enum" types used for whether the basis is
# Real or Complex.  AnyBasis is used when the answer
# is unknown.
##

immutable RealBasis end
immutable ComplexBasis end
immutable AnyBasis end


Base.promote_rule(::Type{RealBasis},::Type{ComplexBasis})=ComplexBasis
Base.promote_rule(::Type{ComplexBasis},::Type{AnyBasis})=AnyBasis
Base.promote_rule(::Type{RealBasis},::Type{AnyBasis})=AnyBasis

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
# D is the domain
# d is the dimension
abstract Space{T,D,d}



typealias RealSpace{D,d} Space{RealBasis,D,d}
typealias ComplexSpace{D,d} Space{ComplexBasis,D,d}
typealias UnivariateSpace{T,D} Space{T,D,1}
typealias BivariateSpace{T} Space{T,AnyDomain,2}
typealias RealUnivariateSpace{D} RealSpace{D,1}




Base.eltype{T}(::Space{T}) = eltype(T)
Base.eltype{T,D,d}(::Type{Space{T,D,d}}) = eltype(T)
basistype{T}(::Space{T}) = T
basistype{T,D,d}(::Type{Space{T,D,d}}) = T
basistype{FT<:Space}(::Type{FT}) = basistype(supertype(FT))

domaintype{T,D}(::Space{T,D}) = D
domaintype{T,D,d}(::Type{Space{T,D,d}}) = D
domaintype{FT<:Space}(::Type{FT}) = domaintype(supertype(FT))

coefficient_type{S}(::Space{S},T) = coefficient_type(S,T)

Base.ndims{S,D,d}(::Space{S,D,d}) = d
dimension(::Space) = ∞  # We assume infinite-dimensional spaces


Space{D<:Number}(d::AbstractVector{D}) = Space(convert(Domain,d))


abstract AmbiguousSpace <: Space{RealBasis,AnyDomain,1}

domain(::AmbiguousSpace) = AnyDomain()


function setdomain{T,D<:Domain}(sp::Space{T,D},d::D)
    S=typeof(sp)
    @assert length(fieldnames(S))==1
    S(d)
end

# function setdomain(sp::Space,d::Domain)
#     S=typeof(sp)
#     @assert length(fieldnames(S))==1
#     # the domain is not compatible, but maybe we c
#     # can drop the space depence.  For example,
#     # CosSpace{Circle{Float64}} -> CosSpace
#     eval(parse(string(S.name.module)*"."*string(S.name)))(d)
# end

setcanonicaldomain(s)=setdomain(s,canonicaldomain(s))
reverseorientation(S::Space)=setdomain(S,reverse(domain(S)))

# AnySpace dictates that an operator can act on any space
# UnsetSpace dictates that an operator is not defined until
#   its domainspace is promoted
# NoSpace is used to indicate no space exists for, e.g.,
# conversion_type

immutable AnySpace <: AmbiguousSpace end
immutable UnsetSpace <: AmbiguousSpace end
immutable NoSpace <: AmbiguousSpace end
immutable ZeroSpace <: AmbiguousSpace end   # ZeroSpace is compatible with all spaces


dimension(::ZeroSpace)=0


isambiguous(::)=false
isambiguous(::AmbiguousSpace) = true

#TODO: should it default to canonicalspace?
points(d::Space,n) = points(domain(d),n)



canonicalspace(T) = T
canonicaldomain(S::Space) = canonicaldomain(domain(S))


# Check whether spaces are the same, override when you need to check parameters
# This is used in place of == to support AnyDomain
spacescompatible{D<:Space}(f::D,g::D) = error("Override spacescompatible for "*string(D))
spacescompatible(::AnySpace,::AnySpace) = true
spacescompatible(::UnsetSpace,::UnsetSpace) = true
spacescompatible(::NoSpace,::NoSpace) = true
spacescompatible(::ZeroSpace,::ZeroSpace) = true
spacescompatible(f,g) = false
==(A::Space,B::Space) = spacescompatible(A,B)&&domain(A)==domain(B)
spacesequal(A::Space,B::Space) = A==B

# check a list of spaces for compatibility
for OP in (:spacescompatible,:domainscompatible,:spacesequal)
    @eval begin
        function $OP{T<:Space}(v::Vector{T})
            for k=1:length(v)-1
                if !$OP(v[k],v[k+1])
                    return false
                end
            end
            true
        end

        $OP{T<:Space}(v::Array{T}) = $OP(vec(v))
    end
end



domain(A::Space) = A.domain # assume it has a field domain



for op in (:tocanonical,:fromcanonical,:tocanonicalD,:fromcanonicalD,:invfromcanonicalD)
    @eval ($op)(sp::Space,x...)=$op(domain(sp),x...)
end

mappoint(a::Space,b::Space,x)=mappoint(domain(a),domain(b),x)
mappoint(a::Space,b::Domain,x)=mappoint(domain(a),b,x)
mappoint(a::Domain,b::Space,x)=mappoint(a,domain(b),x)



for FUNC in (:conversion_rule,:maxspace_rule,:union_rule)
    @eval begin
        function $FUNC(a,b)
            if spacescompatible(a,b)
                a
            else
                NoSpace()
            end
        end
    end
end


for FUNC in (:conversion_type,:maxspace)
    @eval begin
        $FUNC(::AnySpace,::UnsetSpace)=UnsetSpace()
        $FUNC(::UnsetSpace,::AnySpace)=UnsetSpace()
        $FUNC(::ZeroSpace,::UnsetSpace)=UnsetSpace()
        $FUNC(::UnsetSpace,::ZeroSpace)=UnsetSpace()
        $FUNC(::AnySpace,::ZeroSpace)=AnySpace()
        $FUNC(::ZeroSpace,::AnySpace)=AnySpace()
    end

    for TYP in (:AnySpace,:UnsetSpace)
        @eval begin
            $FUNC(a::$TYP,b::$TYP)=a
            $FUNC(a::$TYP,b::Space)=b
            $FUNC(a::Space,b::$TYP)=a
        end
    end
end


# gives a space c that has a banded conversion operator TO a and b
function conversion_type(a,b)
    if spacescompatible(a,b)
        a
    elseif !domainscompatible(a,b)
        NoSpace()  # this avoids having to check eachtime
    else
        cr=conversion_rule(a,b)
        cr==NoSpace()?conversion_rule(b,a):cr
    end
end







# gives a space c that has a banded conversion operator FROM a and b
maxspace(a,b)=NoSpace()  # TODO: this fixes weird bug with Nothing
function maxspace(a::Space,b::Space)
    if spacescompatible(a,b)
        return a
    elseif !domainscompatible(a,b)
        return NoSpace()  # this avoids having to check eachtime
    end



    cr=maxspace_rule(a,b)
    if !isa(cr,NoSpace)
        return cr
    end

    cr=maxspace_rule(b,a)
    if !isa(cr,NoSpace)
        return cr
    end

    cr=conversion_type(a,b)
    if cr==a
        return b
    elseif cr ==b
        return a
    end

    # check if its banded through canonicalspace
    cspa=canonicalspace(a)
    if spacescompatible(cspa,b)
        # we can't call maxspace(cspa,a)
        # maxspace/conversion_type should be implemented for canonicalspace
        error("Override conversion_type or maxspace for "*string(a)*" and "*string(b))
    end
    if cspa != a && maxspace(cspa,a)==cspa
        return maxspace(b,cspa)
    end

    cspb=canonicalspace(b)
    if spacescompatible(cspb,a)
        # we can't call maxspace(cspb,b)
        error("Override conversion_type or maxspace for "*string(a)*" and "*string(b))
    end
    if cspb !=b && maxspace(cspb,b)==cspb
        return maxspace(a,cspb)
    end

    NoSpace()
end




# union combines two spaces
# this is used primarily for addition of two funs
# that may be incompatible
Base.union(a::AmbiguousSpace,b::AmbiguousSpace)=b
Base.union(a::AmbiguousSpace,b::Space)=b
Base.union(a::Space,b::AmbiguousSpace)=a
function Base.union(a::Space,b::Space)
    if spacescompatible(a,b)
        return a
    end

    cr=union_rule(a,b)
    if !isa(cr,NoSpace)
        return cr
    end

    cr=union_rule(b,a)
    if !isa(cr,NoSpace)
        return cr
    end

    cspa=canonicalspace(a)
    cspb=canonicalspace(b)
    if cspa!=a || cspb!=b
        cr=union(cspa,cspb)
    end
    if !isa(cr,NoSpace)
        return cr
    end

    #TODO: Uncomment when Julia bug is fixed
    # cr=maxspace(a,b)  #Max space since we can convert both to it
    # if !isa(cr,NoSpace)
    #     return cr
    # end

    a⊕b
end


# tests whether a Conversion operator exists
hasconversion(a,b)=maxspace(a,b)==b


# tests whether a coefficients can be converted to b
isconvertible(a,b)=hasconversion(a,b)

## Conversion routines
#       coefficients(v::Vector,a,b)
# converts from space a to space b
#       coefficients(v::Fun,a)
# is equivalent to coefficients(v.coefficients,v.space,a)
#       coefficients(v::Vector,a,b,c)
# uses an intermediate space b

coefficients(f,sp1,sp2,sp3)=coefficients(coefficients(f,sp1,sp2),sp2,sp3)

coefficients{T1<:Space,T2<:Space}(f::Vector,::Type{T1},::Type{T2})=coefficients(f,T1(),T2())
coefficients{T1<:Space}(f::Vector,::Type{T1},sp2::Space)=coefficients(f,T1(),sp2)
coefficients{T2<:Space}(f::Vector,sp1::Space,::Type{T2})=coefficients(f,sp1,T2())

## coefficients defaults to calling Conversion, otherwise it tries to pipe through Chebyshev


function defaultcoefficients(f,a,b)
    ct=conversion_type(a,b) # gives a space that has a banded conversion to both a and b

    if spacescompatible(a,b)
        f
    elseif spacescompatible(ct,a)
        (Conversion(a,b)*f).coefficients
    elseif spacescompatible(ct,b)
        (Conversion(b,a)\f).coefficients
    else
        csp=canonicalspace(a)

        if spacescompatible(a,csp)# a is csp, so try b
            csp=canonicalspace(b)
        end
        if spacescompatible(a,csp)||spacescompatible(b,csp)
            # b is csp too, so we are stuck, try Fun constructor
            if domain(b)⊆domain(a)
                coefficients(Fun(x->Fun(f,a)(x),b))
            else
                # we set the value to be zero off the domain of definition
                # but first ensure that domain(b) has a jump
                # TODO: this is disabled as it breaks the case of splitting
                #       one interval into two
                d=domain(a)
#                 if !issubcomponent(d,domain(b))
#                     error("$(d) is not a subcomponent of $(domain(b))")
#                 end

                coefficients(Fun(x->x∈d?Fun(f,a)(x):zero(Fun(f,a)(x)),b))
            end
        else
            coefficients(f,a,csp,b)
        end
    end
end

coefficients(f,a,b)=defaultcoefficients(f,a,b)





## TODO: remove zeros
Base.zero(S::Space)=zeros(S)
Base.zero{T<:Number}(::Type{T},S::Space)=zeros(T,S)
Base.zeros{T<:Number}(::Type{T},S::Space)=Fun(zeros(T,1),S)
Base.zeros(S::Space)=Fun(zeros(1),S)

# catch all
Base.ones(S::Space)=Fun(x->1.0,S)
Base.ones{T<:Number}(::Type{T},S::Space)=Fun(x->one(T),S)

identity_fun(S::Space)=identity_fun(domain(S))

function identity_fun(d::Domain)
    cd=canonicaldomain(d)
    if typeof(d)==typeof(cd)
        Fun(x->x,d) # fall back to constructor
    else
        # this allows support for singularities, that the constructor doesn't
        sf=fromcanonical(d,Fun(identity,cd))
        Fun(coefficients(sf),setdomain(space(sf),d))
    end
end








## rand
# checkpoints is used to give a list of points to double check
# the expansion
Base.rand(d::Space,k...)=rand(domain(d),k...)
checkpoints(d::Space)=checkpoints(domain(d))



## default transforms

# transform converts from values at points(S,n) to coefficients
# itransform converts from coefficients to values at points(S,n)

transform(S::Space,vals)=transform(S,vals,plan_transform(S,vals))
itransform(S::Space,cfs)=itransform(S,cfs,plan_itransform(S,cfs))

function transform(S::Space,vals,plan)
    csp=canonicalspace(S)
    if S==csp
        error("Override transform(::"*string(typeof(S))*",vals,plan)")
    end

    coefficients(transform(csp,vals,plan),csp,S)
end

function itransform(S::Space,cfs,plan)
    csp=canonicalspace(S)
    if S==csp
        error("Override itransform(::"*string(typeof(S))*",cfs,plan)")
    end

    itransform(csp,coefficients(cfs,S,csp),plan)
end


for OP in (:plan_transform,:plan_itransform)
    # plan transform expects a vector
    # this passes an empty Float64 array
    @eval $OP(S::Space,n::Integer)=$OP(S,Array(Float64,n))
end

function plan_transform(S::Space,vals)
    csp=canonicalspace(S)
    if S==csp
        identity #TODO: why identity?
    else
        plan_transform(csp,vals)
    end
end

function plan_itransform(S::Space,cfs)
    csp=canonicalspace(S)
    if S==csp
        identity #TODO: why identity?
    else
        plan_itransform(csp,cfs)
    end
end


## sorting
# we sort spaces lexigraphically by default

for OP in (:<,:(<=),:>,:(>=),:(Base.isless))
    @eval $OP(a::Space,b::Space)=$OP(string(a),string(b))
end
