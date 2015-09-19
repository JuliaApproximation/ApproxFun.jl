immutable ConstantSpace <: UnivariateSpace{RealBasis,AnyDomain} end

ConstantSpace(::AnyDomain)=ConstantSpace()

Fun(c::Number)=Fun([c],ConstantSpace())
Fun(c::Number,d::ConstantSpace)=Fun([c],d)

domain(::ConstantSpace)=AnyDomain()
canonicalspace(C::ConstantSpace)=C
spacescompatible(::ConstantSpace,::ConstantSpace)=true

Base.ones(S::ConstantSpace)=Fun(ones(1),S)
Base.ones(S::Union(AnyDomain,AnySpace,UnsetSpace))=ones(ConstantSpace())
Base.zeros(S::Union(AnyDomain,AnySpace,UnsetSpace))=zeros(ConstantSpace())
evaluate(f::Fun{ConstantSpace},x...)=f.coefficients[1]



promoterangespace(P::Functional,::ConstantSpace,::ConstantSpace)=P # functionals always map to vector space

## Promotion: Zero operators are the only operators that also make sense as functionals
promoterangespace(op::ZeroOperator,::ConstantSpace)=ZeroFunctional(domainspace(op))


# When the union of A and B is a ConstantSpace, then it contains a one
conversion_rule(A::ConstantSpace,B::UnsetSpace)=NoSpace()
conversion_rule(A::ConstantSpace,B::Space)=(union_rule(A,B)==B||union_rule(B,A)==B)?A:NoSpace()


bandinds{S<:Space}(C::Conversion{ConstantSpace,S})=1-length(ones(rangespace(C))),0
function addentries!{S<:Space}(C::Conversion{ConstantSpace,S},A,kr::Range)
    on=ones(rangespace(C))
    for k=kr
        if k≤length(on)
            A[k,1]+=on.coefficients[k]
        end
    end
    A
end

bandinds{F<:Space,T}(D::Multiplication{F,ConstantSpace,T}) = 1-length(D.f),0
function addentries!{F<:Space,T}(D::Multiplication{F,ConstantSpace,T},A,kr)
    Op = Multiplication(D.f,space(D.f))
    for k=kr
        if k≤length(D.f)
            A[k,1]+=Op[k,1]
        end
    end
    A
end
rangespace{F<:Space,T}(D::Multiplication{F,ConstantSpace,T}) = rangespace(Multiplication(D.f,space(D.f)))


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


function addentries!(FO::FunctionalOperator,A,kr::Range)
    if in(1,kr)
        dat=FO.func[1:datalength(FO.func)]
        for j=1:length(dat)
            A[1,j]+=dat[j]
        end
    end
    A
end


for OP in (:+,:-)
    @eval $OP(A::BandedOperator,B::Functional)=$OP(A,FunctionalOperator(B))
    @eval $OP(A::Functional,B::BandedOperator)=$OP(FunctionalOperator(A),B)
end

*(A::BandedOperator,B::Functional)=A*FunctionalOperator(B)
