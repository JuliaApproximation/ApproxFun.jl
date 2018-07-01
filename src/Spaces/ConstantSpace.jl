## Sequence space defintions

# A Fun for SequenceSpace can be an iterator
Base.start(::Fun{SequenceSpace}) = 1
Base.next(f::Fun{SequenceSpace},st) = f[st],st+1
Base.done(f::Fun{SequenceSpace},st) = false # infinite length

getindex(f::Fun{SequenceSpace},k::Integer) =
    k ≤ ncoefficients(f) ? f.coefficients[k] : zero(eltype(f))
getindex(f::Fun{SequenceSpace},K::CartesianIndex{0}) = f[1]
getindex(f::Fun{SequenceSpace},K) = eltype(f)[f[k] for k in K]

Base.length(f::Fun{SequenceSpace}) = ∞


dotu(f::Fun{SequenceSpace},g::Fun{SequenceSpace}) =
    mindotu(f.coefficients,g.coefficients)
dotu(f::Fun{SequenceSpace},g::AbstractVector) =
    mindotu(f.coefficients,g)
dotu(f::AbstractVector,g::Fun{SequenceSpace}) =
    mindotu(f,g.coefficients)

norm(f::Fun{SequenceSpace}) = norm(f.coefficients)
norm(f::Fun{SequenceSpace},k::Int) = norm(f.coefficients,k)
norm(f::Fun{SequenceSpace},k::Number) = norm(f.coefficients,k)


Fun(cfs::AbstractVector,S::SequenceSpace) = Fun(S,cfs)
coefficients(cfs::AbstractVector,::SequenceSpace) = cfs  # all vectors are convertible to SequenceSpace



## Constant space defintions

# setup conversions for spaces that contain constants
macro containsconstants(SP)
    quote
        ApproxFun.union_rule(A::(ApproxFun.ConstantSpace),B::$SP) = B
        Base.promote_rule(A::Type{<:(ApproxFun.ConstantSpace)},B::Type{<:($SP)}) = B

        Base.promote_rule(::Type{ApproxFun.Fun{S,V,VV}},::Type{T}) where {T<:Number,S<:$SP,V,VV} =
            ApproxFun.VFun{S,promote_type(V,T)}
        Base.promote_rule(::Type{ApproxFun.Fun{S}},::Type{T}) where {T<:Number,S<:$SP} = ApproxFun.VFun{S,T}
        Base.promote_rule(::Type{ApproxFun.Fun{S,V,VV}},
        ::Type{Fun{ApproxFun.ConstantSpace{ApproxFun.AnyDomain},T,VT}}) where {T,S<:$SP,V,VV,VT} =
            ApproxFun.VFun{S,promote_type(V,T)}
    end
end



Fun(c::Number) = Fun(ConstantSpace(typeof(c)),[c])
Fun(c::Number,d::ConstantSpace) = Fun(d,[c])

dimension(::ConstantSpace) = 1

#TODO: Change
setdomain(f::Fun{CS},d::Domain) where {CS<:AnyDomain} = Number(f)*one(d)

canonicalspace(C::ConstantSpace) = C
spacescompatible(a::ConstantSpace,b::ConstantSpace)=domainscompatible(a,b)

one(S::ConstantSpace)=Fun(S,fill(1.0,1))
one(S::Union{AnyDomain,UnsetSpace})=one(ConstantSpace())
Base.zeros(S::Union{AnyDomain,UnsetSpace})=zeros(ConstantSpace())
evaluate(f::AbstractVector,::ConstantSpace,x...)=f[1]
evaluate(f::AbstractVector,::ZeroSpace,x...)=zero(eltype(f))


convert(::Type{T}, f::Fun{CS}) where {CS<:ConstantSpace,T<:Number} =
    convert(T, f.coefficients[1])

Number(f::Fun) = convert(Number, f)


# promoting numbers to Fun
# override promote_rule if the space type can represent constants
Base.promote_rule(::Type{Fun{CS}},::Type{T}) where {CS<:ConstantSpace,T<:Number} = Fun{CS,T}
Base.promote_rule(::Type{Fun{CS,V}},::Type{T}) where {CS<:ConstantSpace,T<:Number,V} =
    Fun{CS,promote_type(T,V)}
Base.promote_rule(::Type{IF},::Type{T}) where {T<:Number,IF<:Fun} = Fun


# we know multiplication by constants preserves types
Base.promote_op(::typeof(*),::Type{Fun{CS,T,VT}},::Type{F}) where {CS<:ConstantSpace,T,VT,F<:Fun} =
    promote_op(*,T,F)
Base.promote_op(::typeof(*),::Type{F},::Type{Fun{CS,T,VT}}) where {CS<:ConstantSpace,T,VT,F<:Fun} =
    promote_op(*,F,T)



Base.promote_op(::typeof(LinearAlgebra.matprod),::Type{Fun{S1,T1,VT1}},::Type{Fun{S2,T2,VT2}}) where {S1<:ConstantSpace,T1,VT1,S2<:ConstantSpace,T2,VT2} =
            VFun{promote_type(S1,S2),promote_type(T1,T2)}
Base.promote_op(::typeof(LinearAlgebra.matprod),::Type{Fun{S1,T1,VT1}},::Type{Fun{S2,T2,VT2}}) where {S1<:ConstantSpace,T1,VT1,S2,T2,VT2} =
            VFun{S2,promote_type(T1,T2)}
Base.promote_op(::typeof(LinearAlgebra.matprod),::Type{Fun{S1,T1,VT1}},::Type{Fun{S2,T2,VT2}}) where {S1,T1,VT1,S2<:ConstantSpace,T2,VT2} =
            VFun{S1,promote_type(T1,T2)}


# When the union of A and B is a ConstantSpace, then it contains a one
conversion_rule(A::ConstantSpace,B::UnsetSpace)=NoSpace()
conversion_rule(A::ConstantSpace,B::Space)=(union_rule(A,B)==B||union_rule(B,A)==B) ? A : NoSpace()

conversion_rule(A::ZeroSpace,B::Space) = A
maxspace_rule(A::ZeroSpace,B::Space) = B

Conversion(A::ZeroSpace,B::ZeroSpace) = ConversionWrapper(ZeroOperator(A,B))
Conversion(A::ZeroSpace,B::Space) = ConversionWrapper(ZeroOperator(A,B))

# TODO: this seems like it needs more thought
union_rule(A::ConstantSpace,B::Space) = ConstantSpace(domain(B))⊕B


## Special Multiplication and Conversion for constantspace

#  TODO: this is a special work around but really we want it to be blocks
Conversion(a::ConstantSpace,b::Space{D}) where {D<:BivariateDomain} = ConcreteConversion{typeof(a),typeof(b),
        promote_type(real(prectype(a)),real(prectype(b)))}(a,b)

Conversion(a::ConstantSpace,b::Space) = ConcreteConversion(a,b)
bandinds(C::ConcreteConversion{CS,S}) where {CS<:ConstantSpace,S<:Space} =
    1-ncoefficients(one(rangespace(C))),0
function getindex(C::ConcreteConversion{CS,S,T},k::Integer,j::Integer) where {CS<:ConstantSpace,S<:Space,T}
    if j != 1
        throw(BoundsError())
    end
    on=one(rangespace(C))
    k ≤ ncoefficients(on) ? convert(T,on.coefficients[k]) : zero(T)
end

coefficients(f::AbstractVector,sp::ConstantSpace{Segment{Vec{2,TT}}},
             ts::TensorSpace{SV,DD}) where {TT,SV,DD<:BivariateDomain} =
    f[1]*one(ts).coefficients
coefficients(f::AbstractVector,sp::ConstantSpace,ts::Space) = f[1]*one(ts).coefficients


########
# Evaluation
########

#########
# Multiplication
#########


# this is identity operator, but we don't use MultiplicationWrapper to avoid
# ambiguity errors

defaultMultiplication(f::Fun{CS},b::ConstantSpace) where {CS<:ConstantSpace} =
    ConcreteMultiplication(f,b)
defaultMultiplication(f::Fun{CS},b::Space) where {CS<:ConstantSpace} =
    ConcreteMultiplication(f,b)
defaultMultiplication(f::Fun,b::ConstantSpace) = ConcreteMultiplication(f,b)

bandinds(D::ConcreteMultiplication{CS1,CS2,T}) where {CS1<:ConstantSpace,CS2<:ConstantSpace,T} =
    0,0
getindex(D::ConcreteMultiplication{CS1,CS2,T},k::Integer,j::Integer) where {CS1<:ConstantSpace,CS2<:ConstantSpace,T} =
    k==j==1 ? convert(T,D.f.coefficients[1]) : one(T)

rangespace(D::ConcreteMultiplication{CS1,CS2,T}) where {CS1<:ConstantSpace,CS2<:ConstantSpace,T} =
    D.space


rangespace(D::ConcreteMultiplication{F,UnsetSpace,T}) where {F<:ConstantSpace,T} =
    UnsetSpace()
bandinds(D::ConcreteMultiplication{F,UnsetSpace,T}) where {F<:ConstantSpace,T} =
    (-∞,∞)
getindex(D::ConcreteMultiplication{F,UnsetSpace,T},k::Integer,j::Integer) where {F<:ConstantSpace,T} =
    error("No range space attached to Multiplication")



bandinds(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = 0,0
blockbandinds(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = 0,0
subblockbandinds(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = 0,0
subblockbandinds(D::ConcreteMultiplication{CS,F,T}, k) where {CS<:ConstantSpace,F<:Space,T} = 0
isbandedblockbanded(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = true
isblockbanded(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = true
getindex(D::ConcreteMultiplication{CS,F,T},k::Integer,j::Integer) where {CS<:ConstantSpace,F<:Space,T} =
    k==j ? convert(T, D.f) : zero(T)
rangespace(D::ConcreteMultiplication{CS,F,T}) where {CS<:ConstantSpace,F<:Space,T} = D.space


bandinds(D::ConcreteMultiplication{F,CS,T}) where {CS<:ConstantSpace,F<:Space,T} = 1-ncoefficients(D.f),0
function getindex(D::ConcreteMultiplication{F,CS,T},k::Integer,j::Integer) where {CS<:ConstantSpace,F<:Space,T}
    k≤ncoefficients(D.f) && j==1 ? convert(T,D.f.coefficients[k]) : zero(T)
end
rangespace(D::ConcreteMultiplication{F,CS,T}) where {CS<:ConstantSpace,F<:Space,T} = space(D.f)



# functionals always map to Constant space
function promoterangespace(P::Operator,A::ConstantSpace,cur::ConstantSpace)
    @assert isafunctional(P)
    domain(A)==domain(cur) ? P : SpaceOperator(P,domainspace(P),A)
end


for op = (:*,:/)
    @eval $op(f::Fun,c::Fun{CS}) where {CS<:ConstantSpace} = f*convert(Number,c)
end



## Multivariate case
union_rule(a::TensorSpace,b::ConstantSpace{AnyDomain})=TensorSpace(map(sp->union(sp,b),a.spaces))
## Special spaces

function convert(::Type{T},f::Fun{TS}) where {TS<:TensorSpace,T<:Number}
    if all(sp->isa(sp,ConstantSpace),space(f).spaces)
        convert(T,f.coefficients[1])
    else
        error("Cannot convert $f to type $T")
    end
end

convert(::Type{T},
            f::Fun{TensorSpace{Tuple{CS1,CS2},DD,RR}}) where {CS1<:ConstantSpace,CS2<:ConstantSpace,T<:Number,DD,RR} =
    convert(T,f.coefficients[1])

isconstspace(sp::TensorSpace) = all(isconstspace,sp.spaces)


# Supports constants in operators
promoterangespace(M::ConcreteMultiplication{CS,UnsetSpace},
                             ps::UnsetSpace) where {CS<:ConstantSpace} = M
promoterangespace(M::ConcreteMultiplication{CS,UnsetSpace},
                             ps::Space) where {CS<:ConstantSpace} =
                        promoterangespace(Multiplication(M.f,space(M.f)),ps)

# Possible hack: we try uing constant space for [1 Operator()] \ z.
choosedomainspace(M::ConcreteMultiplication{D,UnsetSpace},sp::UnsetSpace) where {D<:ConstantSpace} = space(M.f)
choosedomainspace(M::ConcreteMultiplication{D,UnsetSpace},sp::Space) where {D<:ConstantSpace} = space(M.f)

Base.isfinite(f::Fun{CS}) where {CS<:ConstantSpace} = isfinite(Number(f))
