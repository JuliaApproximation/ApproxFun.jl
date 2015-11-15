

"""
`ConstantSpace` Represents a single number.  The remaining
coefficients are ignored.
"""

immutable ConstantSpace <: UnivariateSpace{RealBasis,AnyDomain} end

ConstantSpace(::AnyDomain)=ConstantSpace()

Fun(c::Number)=Fun([c],ConstantSpace())
Fun(c::Number,d::ConstantSpace)=Fun([c],d)

dimension(::ConstantSpace)=1

domain(::ConstantSpace)=AnyDomain()

setdomain(f::Fun{ConstantSpace},d::Domain)=Number(f)*ones(d)

canonicalspace(C::ConstantSpace)=C
spacescompatible(::ConstantSpace,::ConstantSpace)=true

Base.ones(S::ConstantSpace)=Fun(ones(1),S)
Base.ones(S::Union{AnyDomain,AnySpace,UnsetSpace})=ones(ConstantSpace())
Base.zeros(S::Union{AnyDomain,AnySpace,UnsetSpace})=zeros(ConstantSpace())
evaluate(f::AbstractVector,::ConstantSpace,x...)=f[1]
evaluate(f::AbstractVector,::ConstantSpace,x::Array)=f[1]*ones(x)

evaluate(f::AbstractVector,::ZeroSpace,x...)=zero(eltype(f))
evaluate(f::AbstractVector,::ZeroSpace,x::Array)=zeros(x)


Base.convert{T<:Number}(::Type{T},f::Fun{ConstantSpace})=convert(T,f.coefficients[1])

# promoting numbers to Fun
# override promote_rule if the space type can represent constants
Base.promote_rule{T<:Number}(::Type{Fun{ConstantSpace}},::Type{T})=Fun{ConstantSpace,T}
Base.promote_rule{T<:Number,V}(::Type{Fun{ConstantSpace,V}},::Type{T})=Fun{ConstantSpace,promote_type(T,V)}
Base.promote_rule{T<:Number,IF<:Fun}(::Type{IF},::Type{T})=Fun



promoterangespace(P::Functional,::ConstantSpace,::ConstantSpace)=P # functionals always map to vector space

## Promotion: Zero operators are the only operators that also make sense as functionals
promoterangespace(op::ZeroOperator,::ConstantSpace)=ZeroFunctional(domainspace(op))


# When the union of A and B is a ConstantSpace, then it contains a one
conversion_rule(A::ConstantSpace,B::UnsetSpace)=NoSpace()
conversion_rule(A::ConstantSpace,B::Space)=(union_rule(A,B)==B||union_rule(B,A)==B)?A:NoSpace()

conversion_rule(A::ZeroSpace,B::Space)=A
maxspace_rule(A::ZeroSpace,B::Space)=B
Conversion(A::ZeroSpace,B::Space)=ConversionWrapper(SpaceOperator(ZeroOperator(),A,B))


bandinds{S<:Space}(C::Conversion{ConstantSpace,S})=1-length(ones(rangespace(C))),0
function addentries!{S<:Space}(C::Conversion{ConstantSpace,S},A,kr::Range,::Colon)
    on=ones(rangespace(C))
    for k=kr
        if k≤length(on)
            A[k,1]+=on.coefficients[k]
        end
    end
    A
end

bandinds{F<:Space,T}(D::Multiplication{F,ConstantSpace,T}) = 1-length(D.f),0
function addentries!{F<:Space,T}(D::Multiplication{F,ConstantSpace,T},A,kr::Range,::Colon)
    Op = Multiplication(D.f,space(D.f))
    for k=kr
        if k≤length(D.f)
            A[k,1]+=Op[k,1]
        end
    end
    A
end

function addentries!{T}(D::Multiplication{ConstantSpace,ConstantSpace,T},A,kr::Range,::Colon)
    if 1 in kr
        A[1,1]+=D.f.coefficients[1]
    end
    A
end

rangespace{F<:Space,T}(D::Multiplication{F,ConstantSpace,T}) = rangespace(Multiplication(D.f,space(D.f)))
rangespace{T}(D::Multiplication{ConstantSpace,ConstantSpace,T}) = ConstantSpaec()


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

*{T,D<:Union{DefiniteIntegral,DefiniteLineIntegral},M<:AbstractMultiplication,V}(A::FunctionalOperator{TimesFunctional{T,D,M},V},b::Fun) = Fun(A.func*b)
