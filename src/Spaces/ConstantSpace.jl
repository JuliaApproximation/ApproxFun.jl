

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

dimension(::ConstantSpace)=1

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


# functionals always map to Constant space
promoterangespace(P::Functional,A::ConstantSpace,cur::ConstantSpace)=domain(A)==domain(cur)?P:SpaceFunctional(P,domainspace(P),domain(A))

## Promotion: Zero operators are the only operators that also make sense as functionals
promoterangespace(op::ZeroOperator,::ConstantSpace)=ZeroFunctional(domainspace(op))


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
bandinds{CS<:ConstantSpace,S<:Space}(C::ConcreteConversion{CS,S})=1-length(ones(rangespace(C))),0
function addentries!{CS<:ConstantSpace,S<:Space}(C::ConcreteConversion{CS,S},A,kr::Range,::Colon)
    on=ones(rangespace(C))
    for k=kr
        if k≤length(on)
            A[k,1]+=on.coefficients[k]
        end
    end
    A
end


# this is identity operator, but we don't use MultiplicationWrapper to avoid
# ambiguity errors


bandinds{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T}) = 0,0
function addentries!{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T},A,kr::Range,::Colon)
    if 1 in kr
        A[1,1]+=D.f.coefficients[1]
    end
    A
end
rangespace{CS1<:ConstantSpace,CS2<:ConstantSpace,T}(D::ConcreteMultiplication{CS1,CS2,T}) = D.space


rangespace{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T})=UnsetSpace()
bandinds{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T})=error("No range space attached to Multiplication")
addentries!{F<:ConstantSpace,T}(D::ConcreteMultiplication{F,UnsetSpace,T},A,kr::Range,::Colon)=error("No range space attached to Multiplication")



bandinds{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T}) = 0,0
function addentries!{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T},A,kr::Range,::Colon)
    c=Number(D.f)
    for k=kr
        A[k,k]+=c
    end
    A
end
rangespace{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{CS,F,T}) = D.space


bandinds{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T}) = 1-length(D.f),0
function addentries!{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T},A,kr::Range,::Colon)
    Op = Multiplication(D.f,space(D.f))
    for k=kr
        if k≤length(D.f)
            A[k,1]+=Op[k,1]
        end
    end
    A
end
rangespace{CS<:ConstantSpace,F<:Space,T}(D::ConcreteMultiplication{F,CS,T}) = rangespace(Multiplication(D.f,space(D.f)))




###
# FunctionalOperator treats a functional like an operator
###

immutable FunctionalOperator{FT,T} <: BandedOperator{T}
    func::FT
end

FunctionalOperator(func::Functional)=FunctionalOperator{typeof(func),eltype(func)}(func)

function Base.convert{O<:Operator}(::Type{O},FO::FunctionalOperator)
    if eltype(O)==eltype(FO)
        FO::O
    else
        FunctionalOperator{typeof(FO.func),eltype(O)}(FO.func)
    end
end

bandinds(FO::FunctionalOperator)=0,datalength(FO.func)-1
domainspace(FO::FunctionalOperator)=domainspace(FO.func)
rangespace(FO::FunctionalOperator)=ConstantSpace()

for TYP in (:AnySpace,:UnsetSpace,:Space)
    @eval promotedomainspace(FT::FunctionalOperator,sp::$TYP)=FunctionalOperator(promotedomainspace(FT.func,sp))
end


function addentries!(FO::FunctionalOperator,A,kr::Range,::Colon)
    if in(1,kr)
        dat=FO.func[1:datalength(FO.func)]
        for j=1:length(dat)
            A[1,j]+=dat[j]
        end
    end
    A
end

function *(f::Fun,A::Functional)
    if datalength(A)<Inf
        # We get a BandedOperator, so we take that into account
        TimesOperator(Multiplication(f,ConstantSpace()),FunctionalOperator(A))
    else
        LowRankOperator(f,A)
    end
end

Base.convert(::Type{BandedOperator},B::Functional)=FunctionalOperator(B)
Base.convert(::Type{BandedBelowOperator},B::Functional)=datalength(B)<Inf?convert(BandedOperator,B):Fun(1,ConstantSpace())*B



for OP in (:+,:-)
    @eval $OP(A::BandedOperator,B::Functional)=$OP(A,convert(BandedBelowOperator,B))
    @eval $OP(A::Functional,B::BandedOperator)=$OP(convert(BandedBelowOperator,A),B)
end

*(A::BandedOperator,B::Functional)=A*convert(BandedBelowOperator,B)

*{T,D<:Union{DefiniteIntegral,DefiniteLineIntegral},M<:Multiplication,V}(A::FunctionalOperator{TimesFunctional{T,D,M},V},b::Fun) = Fun(A.func*b)


for op = (:*,:.*,:./,:/)
    @eval $op{CS<:ConstantSpace}(f::Fun,c::Fun{CS}) = f*convert(Number,c)
end



## Multivariate case
union_rule(a::TensorSpace,b::ConstantSpace{AnyDomain})=TensorSpace(map(sp->union(sp,b),a.spaces))

