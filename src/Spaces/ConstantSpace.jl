

"""
`ConstantSpace` Represents a single number.  The remaining
coefficients are ignored.
"""

immutable ConstantSpace{DD} <: UnivariateSpace{RealBasis,DD}
    domain::DD
    ConstantSpace(d::DD)=new(d)
end

ConstantSpace(d::Domain)=ConstantSpace{typeof(d)}(d)
ConstantSpace()=ConstantSpace(AnyDomain())

# we override maxspace instead of maxspace_rule to avoid
# domainscompatible check.
for OP in (:maxspace,:(Base.union))
    @eval begin
        $OP(A::ConstantSpace{AnyDomain},B::ConstantSpace{AnyDomain})=A
        $OP(A::ConstantSpace{AnyDomain},B::ConstantSpace)=B
        $OP(A::ConstantSpace,B::ConstantSpace{AnyDomain})=A
        $OP(A::ConstantSpace,B::ConstantSpace)=ConstantSpace(domain(A) ∪ domain(B))
    end
end

Fun(c::Number)=Fun([c],ConstantSpace())
Fun(c::Number,d::ConstantSpace)=Fun([c],d)

dimension(::ConstantSpace) = 1

#TODO: Change
setdomain{CS<:AnyDomain}(f::Fun{CS},d::Domain)=Number(f)*ones(d)

canonicalspace(C::ConstantSpace)=C
spacescompatible(a::ConstantSpace,b::ConstantSpace)=domainscompatible(a,b)

Base.ones(S::ConstantSpace)=Fun(ones(1),S)
Base.ones(S::Union{AnyDomain,AnySpace,UnsetSpace})=ones(ConstantSpace())
Base.zeros(S::Union{AnyDomain,AnySpace,UnsetSpace})=zeros(ConstantSpace())
evaluate(f::AbstractVector,::ConstantSpace,x...)=f[1]
evaluate(f::AbstractVector,::ConstantSpace,x::Array)=f[1]*ones(x)

evaluate(f::AbstractVector,::ZeroSpace,x...)=zero(eltype(f))
evaluate(f::AbstractVector,::ZeroSpace,x::Array)=zeros(x)


Base.convert{CS<:ConstantSpace,T<:Number}(::Type{T},f::Fun{CS})=convert(T,f.coefficients[1])

# promoting numbers to Fun
# override promote_rule if the space type can represent constants
Base.promote_rule{CS<:ConstantSpace,T<:Number}(::Type{Fun{CS}},::Type{T})=Fun{CS,T}
Base.promote_rule{CS<:ConstantSpace,T<:Number,V}(::Type{Fun{CS,V}},::Type{T})=Fun{CS,promote_type(T,V)}
Base.promote_rule{T<:Number,IF<:Fun}(::Type{IF},::Type{T})=Fun


# When the union of A and B is a ConstantSpace, then it contains a one
conversion_rule(A::ConstantSpace,B::UnsetSpace)=NoSpace()
conversion_rule(A::ConstantSpace,B::Space)=(union_rule(A,B)==B||union_rule(B,A)==B)?A:NoSpace()

conversion_rule(A::ZeroSpace,B::Space)=A
maxspace_rule(A::ZeroSpace,B::Space)=B
Conversion(A::ZeroSpace,B::Space)=ConversionWrapper(SpaceOperator(ZeroOperator(),A,B))


union_rule(A::ConstantSpace,B::Space)=ConstantSpace(domain(B))⊕B


## Special Multiplication and Conversion for constantspace

#  TODO: this is a special work around but really we want it to be blocks
Conversion{T,D}(a::ConstantSpace,b::Space{T,D,2})=ConcreteConversion{typeof(a),typeof(b),
        promote_type(op_eltype_realdomain(a),eltype(op_eltype_realdomain(b)))}(a,b)

Conversion(a::ConstantSpace,b::Space)=ConcreteConversion(a,b)
bandinds{CS<:ConstantSpace,S<:Space}(C::ConcreteConversion{CS,S})=1-ncoefficients(ones(rangespace(C))),0
function getindex{CS<:ConstantSpace,S<:Space,T}(C::ConcreteConversion{CS,S,T},k::Integer,j::Integer)
    if j != 1
        throw(BoundsError())
    end
    on=ones(rangespace(C))
    k ≤ ncoefficients(on)?T(on.coefficients[k]):zero(T)
end


# this is identity operator, but we don't use MultiplicationWrapper to avoid
# ambiguity errors

defaultMultiplication{CS<:ConstantSpace}(f::Fun{CS},b::ConstantSpace) =
    ConcreteMultiplication(f,b)
defaultMultiplication{CS<:ConstantSpace}(f::Fun{CS},b::Space) =
    ConcreteMultiplication(f,b)
defaultMultiplication(f::Fun,b::ConstantSpace) = ConcreteMultiplication(f,b)

bandinds{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T}) =
    0,0
getindex{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T},k::Integer,j::Integer) =
    k==j==1?T(D.f.coefficients[1]):one(T)

rangespace{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T}) =
    D.space


rangespace{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T}) =
    UnsetSpace()
bandinds{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T}) =
    error("No range space attached to Multiplication")
getindex{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T},k::Integer,j::Integer) =
    error("No range space attached to Multiplication")



bandinds{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T}) = 0,0
getindex{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T},k::Integer,j::Integer) =
    k==j?T(D.f):zero(T)
rangespace{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T}) = D.space


bandinds{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T}) = 1-ncoefficients(D.f),0
function getindex{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T},k::Integer,j::Integer)
    Op = Multiplication(D.f,space(D.f))
    k≤ncoefficients(D.f) && j==1?T(Op[k,1]):zero(T)
end
rangespace{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T}) =
    rangespace(Multiplication(D.f,space(D.f)))



# functionals always map to Constant space
function promoterangespace(P::Operator,A::ConstantSpace,cur::ConstantSpace)
    @assert isafunctional(P)
    domain(A)==domain(cur)?P:SpaceOperator(P,domainspace(P),A)
end



for op = (:*,:.*,:./,:/)
    @eval $op{CS<:ConstantSpace}(f::Fun,c::Fun{CS}) = f*convert(Number,c)
end



## Multivariate case
union_rule(a::TensorSpace,b::ConstantSpace{AnyDomain})=TensorSpace(map(sp->union(sp,b),a.spaces))
