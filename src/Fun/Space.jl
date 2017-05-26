

export Space, domainspace, rangespace, maxspace,Space,conversion_type, transform,
            itransform, SequenceSpace, ConstantSpace



# Space maps the Domain to the type R
# For example, we have
#   Chebyshev{Interval{Float64}} <: Space{Interval{Float64},Float64}
#   Laurent{PeriodicInterval{Float64}} <: Space{PeriodicInterval{Float64},Complex128}
#   Fourier{Circle{Complex128}} <: Space{Circle{Complex128},Float64}
# Note for now Space doesn't contain any information about the coefficients

abstract type Space{D,R} end



const RealSpace = Space{D,R} where {D,R<:Real}
const ComplexSpace = Space{D,R} where {D,R<:Complex}
const UnivariateSpace = Space{D,R} where {D<:UnivariateDomain,R}
const BivariateSpace = Space{D,R}  where {D<:BivariateDomain,R}
const RealUnivariateSpace = RealSpace{D,R} where {D<:UnivariateDomain,R<:Real}




Base.eltype{T}(::Space{T}) = error("Don't use")
Base.eltype{D,R}(::Type{Space{D,R}}) = error("Don't use")

domaintype{D,R}(::Space{D,R}) = D
domaintype{D,R}(::Type{Space{D,R}}) = D
domaintype{FT<:Space}(::Type{FT}) = domaintype(supertype(FT))
rangetype{D,R}(::Space{D,R}) = R
rangetype{D,R}(::Type{Space{D,R}}) = R
rangetype{FT<:Space}(::Type{FT}) = rangetype(supertype(FT))

domaindimension(sp::Space) = dimension(domain(sp))
dimension(::Space) = ∞  # We assume infinite-dimensional spaces


# add indexing for all spaces, not just DirectSumSpace
# mimicking scalar vs vector

# TODO: 0.5 iterator
Base.start(s::Space) = false
Base.next(s::Space,st) = (s,true)
Base.done(s::Space,st) = st
Base.length(s::Space) = 1
Base.getindex(s::Space,::CartesianIndex{0}) = s
getindex(s::Space,k) = k == 1 ? s : throw(BoundsError())
Base.endof(s::Space) = 1


#supports broadcasting, overloaded for ArraySpace
Base.size(::Space) = ()

# but broadcasts like Number (even for overloaded ArraySpace)
# TODO for v0.6: remove this
Base.broadcast(f,x::Union{Number,Domain,Space}...) = f(x...)

# the default is all spaces have one-coefficient blocks
blocklengths(S::Space) = repeated(true,dimension(S))
nblocks(S::Space) = length(blocklengths(S))
block(S::Space,k) = Block(k)

Space{D<:Number}(d::AbstractVector{D}) = Space(convert(Domain,d))


abstract type AmbiguousSpace <: Space{AnyDomain,UnsetNumber} end
domain(::AmbiguousSpace) = AnyDomain()


function setdomain{D<:Domain}(sp::Space{D},d::D)
    S = typeof(sp)
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

setcanonicaldomain(s) = setdomain(s,canonicaldomain(s))
reverseorientation(S::Space) = setdomain(S,reverse(domain(S)))


# UnsetSpace dictates that an operator is not defined until
#   its domainspace is promoted
# NoSpace is used to indicate no space exists for, e.g.,
# conversion_type

struct UnsetSpace <: AmbiguousSpace end
struct NoSpace <: AmbiguousSpace end
struct ZeroSpace <: AmbiguousSpace end   # ZeroSpace is compatible with all spaces


dimension(::ZeroSpace) = 0


isambiguous(::)=false
isambiguous(::AmbiguousSpace) = true

#TODO: should it default to canonicalspace?
points(d::Space,n) = points(domain(d),n)



canonicalspace(T) = T
canonicaldomain(S::Space) = canonicaldomain(domain(S))


# Check whether spaces are the same, override when you need to check parameters
# This is used in place of == to support AnyDomain
spacescompatible{D<:Space}(f::D,g::D) = error("Override spacescompatible for "*string(D))
spacescompatible(::UnsetSpace,::UnsetSpace) = true
spacescompatible(::NoSpace,::NoSpace) = true
spacescompatible(::ZeroSpace,::ZeroSpace) = true
spacescompatible(f,g) = false
==(A::Space,B::Space) = spacescompatible(A,B)&&domain(A)==domain(B)
spacesequal(A::Space,B::Space) = A==B

# check a list of spaces for compatibility
for OP in (:spacescompatible,:domainscompatible,:spacesequal),TYP in (:AbstractArray,:Tuple)
    @eval function $OP(v::$TYP)
        for k=1:length(v)-1
            if !$OP(v[k],v[k+1])
                return false
            end
        end
        true
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
        $FUNC(::ZeroSpace,::UnsetSpace) = UnsetSpace()
        $FUNC(::UnsetSpace,::ZeroSpace) = UnsetSpace()
        $FUNC(a::UnsetSpace,b::UnsetSpace) = a
        $FUNC(a::UnsetSpace,b::Space) = b
        $FUNC(a::Space,b::UnsetSpace) = a
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
maxspace(a,b) = NoSpace()  # TODO: this fixes weird bug with Nothing
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
        if isambiguous(domain(a))
            return b
        else
            return a
        end
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

Base.union(a::Space) = a

Base.union(a::Space,b::Space,c::Space) = union(union(a,b),c)
Base.union(a::Space,b::Space,c::Space,d::Space...) =
    union(union(a,b),c,d...)


# tests whether a Conversion operator exists
hasconversion(a,b) = maxspace(a,b)==b


# tests whether a coefficients can be converted to b
isconvertible(a,b) = hasconversion(a,b)

## Conversion routines
#       coefficients(v::AbstractVector,a,b)
# converts from space a to space b
#       coefficients(v::Fun,a)
# is equivalent to coefficients(v.coefficients,v.space,a)
#       coefficients(v::AbstractVector,a,b,c)
# uses an intermediate space b

coefficients(f,sp1,sp2,sp3) = coefficients(coefficients(f,sp1,sp2),sp2,sp3)

coefficients{T1<:Space,T2<:Space}(f::AbstractVector,::Type{T1},::Type{T2}) =
    coefficients(f,T1(),T2())
coefficients{T1<:Space}(f::AbstractVector,::Type{T1},sp2::Space) = coefficients(f,T1(),sp2)
coefficients{T2<:Space}(f::AbstractVector,sp1::Space,::Type{T2}) = coefficients(f,sp1,T2())

## coefficients defaults to calling Conversion, otherwise it tries to pipe through Chebyshev


function defaultcoefficients(f,a,b)
    ct=conversion_type(a,b) # gives a space that has a banded conversion to both a and b

    if spacescompatible(a,b)
        f
    elseif spacescompatible(ct,a)
        A_mul_B_coefficients(Conversion(a,b),f)
    elseif spacescompatible(ct,b)
        A_ldiv_B_coefficients(Conversion(b,a),f)
    else
        csp=canonicalspace(a)

        if spacescompatible(a,csp)# a is csp, so try b
            csp=canonicalspace(b)
        end
        if spacescompatible(a,csp)||spacescompatible(b,csp)
            # b is csp too, so we are stuck, try Fun constructor
            if domain(b)⊆domain(a)
                coefficients(Fun(x->Fun(a,f)(x),b))
            else
                # we set the value to be zero off the domain of definition
                # but first ensure that domain(b) has a jump
                # TODO: this is disabled as it breaks the case of splitting
                #       one interval into two
                d=domain(a)
#                 if !issubcomponent(d,domain(b))
#                     error("$(d) is not a subcomponent of $(domain(b))")
#                 end

                coefficients(Fun(x->x∈d?Fun(a,f)(x):zero(Fun(a,f)(x)),b))
            end
        else
            coefficients(f,a,csp,b)
        end
    end
end

coefficients(f,a,b) = defaultcoefficients(f,a,b)





## TODO: remove zeros
Base.zero(S::Space) = zeros(S)
Base.zero{T<:Number}(::Type{T},S::Space) = zeros(T,S)
Base.zeros{T<:Number}(::Type{T},S::Space) = Fun(S,zeros(T,1))
Base.zeros(S::Space) = Fun(S,zeros(1))

# catch all
Base.ones(S::Space) = Fun(x->1.0,S)
Base.ones{T<:Number}(::Type{T},S::Space) = Fun(x->one(T),S)

identity_fun(S::Space) = identity_fun(domain(S))

function identity_fun(d::Domain)
    cd=canonicaldomain(d)
    if typeof(d) == typeof(cd)
        Fun(x->x,d) # fall back to constructor
    else
        # this allows support for singularities, that the constructor doesn't
        sf=fromcanonical(d,Fun(identity,cd))
        Fun(setdomain(space(sf),d),coefficients(sf))
    end
end








## rand
# checkpoints is used to give a list of points to double check
# the expansion
Base.rand(d::Space,k...) = rand(domain(d),k...)
checkpoints(d::Space) = checkpoints(domain(d))



## default transforms

# These plans are use to wrap another plan
for Typ in (:TransformPlan,:ITransformPlan)
    @eval begin
        struct $Typ{T,SP,inplace,PL} <: FFTW.Plan{T}
            space::SP
            plan::PL
        end
        $Typ{inplace}(space,plan,::Type{Val{inplace}}) =
            $Typ{eltype(plan),typeof(space),inplace,typeof(plan)}(space,plan)
    end
end

for Typ in (:CanonicalTransformPlan,:ICanonicalTransformPlan)
    @eval begin
        struct $Typ{T,SP,PL,CSP} <: FFTW.Plan{T}
            space::SP
            plan::PL
            canonicalspace::CSP
        end
        $Typ(space,plan,csp) =
            $Typ{eltype(plan),typeof(space),typeof(plan),typeof(csp)}(space,plan,csp)
        $Typ{T}(::Type{T},space,plan,csp) =
            $Typ{T,typeof(space),typeof(plan),typeof(csp)}(space,plan,csp)
    end
end


# Canonical plan uses coefficients
function CanonicalTransformPlan(space,v)
    csp = canonicalspace(space)
    CanonicalTransformPlan(eltype(v),space,plan_transform(csp,v),csp)
end
function plan_transform(sp::Space,vals)
    csp = canonicalspace(sp)
    if sp == csp
        error("Override for $sp")
    end
    CanonicalTransformPlan(sp,plan_transform(csp,vals),csp)
end

function ICanonicalTransformPlan(space,v)
    csp = canonicalspace(space)
    cfs = coefficients(v,space,csp)
    ICanonicalTransformPlan(eltype(v),space,plan_itransform(csp,cfs),csp)
end
function plan_itransform(sp::Space,v)
    csp = canonicalspace(sp)
    if sp == csp
        error("Override for $sp")
    end
    cfs = coefficients(v,sp,csp)
    ICanonicalTransformPlan(sp,plan_itransform(csp,cfs),csp)
end


plan_transform!(sp::Space,vals) = error("Override for $sp")
plan_itransform!(sp::Space,cfs) = error("Override for $sp")



# transform converts from values at points(S,n) to coefficients
# itransform converts from coefficients to values at points(S,n)

transform(S::Space,vals) = plan_transform(S,vals)*vals
itransform(S::Space,cfs) = plan_itransform(S,cfs)*cfs

*(P::CanonicalTransformPlan,vals::AbstractVector) = coefficients(P.plan*vals,P.canonicalspace,P.space)
*(P::ICanonicalTransformPlan,cfs::AbstractVector) = P.plan*coefficients(cfs,P.space,P.canonicalspace)



for OP in (:plan_transform,:plan_itransform,:plan_transform!,:plan_itransform!)
    # plan transform expects a vector
    # this passes an empty Float64 array
    @eval $OP(S::Space,n::Integer) = $OP(S,Vector{Float64}(n))
end

## sorting
# we sort spaces lexigraphically by default

for OP in (:<,:(<=),:>,:(>=),:(Base.isless))
    @eval $OP(a::Space,b::Space)=$OP(string(a),string(b))
end



## Important special spaces


"""
`ConstantSpace` is the 1-dimensional scalar space.
"""

struct ConstantSpace{DD,R} <: Space{DD,R}
    domain::DD
    ConstantSpace{DD,R}(d::DD) where {DD,R} = new(d)
    ConstantSpace{DD,R}(d::AnyDomain) where {DD,R} = new(DD(d))
end

ConstantSpace(d::Domain) = ConstantSpace{typeof(d),real(prectype(d))}(d)
ConstantSpace() = ConstantSpace(AnyDomain())

isconstspace(::ConstantSpace) = true

for pl in (:plan_transform,:plan_transform!,:plan_itransform,:plan_itransform!)
    @eval $pl(sp::ConstantSpace,vals::AbstractVector) = I
end

# we override maxspace instead of maxspace_rule to avoid
# domainscompatible check.
for OP in (:maxspace,:(Base.union))
    @eval begin
        $OP(A::ConstantSpace{AnyDomain},B::ConstantSpace{AnyDomain}) = A
        $OP(A::ConstantSpace{AnyDomain},B::ConstantSpace) = B
        $OP(A::ConstantSpace,B::ConstantSpace{AnyDomain}) = A
        $OP(A::ConstantSpace,B::ConstantSpace) = ConstantSpace(domain(A) ∪ domain(B))
    end
end



# Range type is Void since function evaluation is not defined
struct SequenceSpace <: Space{PositiveIntegers,Void} end

doc"""
`SequenceSpace` is the space of all sequences, i.e., infinite vectors.
Also denoted ℓ⁰.
"""
SequenceSpace()


const ℓ⁰ = SequenceSpace()
dimension(::SequenceSpace) = ∞
domain(::SequenceSpace) = ℕ
spacescompatible(::SequenceSpace,::SequenceSpace) = true


## Boundary

∂(S::Space) = ∂(domain(S))
