immutable ConstantSpace <: UnivariateSpace{RealBasis} end
domain(::ConstantSpace)=AnyDomain()
canonicalspace(C::ConstantSpace)=C
spacescompatible(::ConstantSpace,::ConstantSpace)=true

promoterangespace(P::Functional,::ConstantSpace,::ConstantSpace)=P # functionals always map to vector space

## Promotion: Zero operators are the only operators that also make sense as functionals
promoterangespace(op::ZeroOperator,::ConstantSpace)=ZeroFunctional(domainspace(op))


# We assume we can always convert a constant to the functionspace
# this is only true if 1 is an element
conversion_rule(A::ConstantSpace,B::FunctionSpace)=A

bandinds{S<:FunctionSpace}(C::Conversion{ConstantSpace,S})=1-length(ones(rangespace(C))),0
function addentries!{S<:FunctionSpace}(C::Conversion{ConstantSpace,S},A,kr::Range)
    on=ones(rangespace(C))
    for k=kr
        if k≤length(on)
            A[k,1]+=on.coefficients[k]
        end
    end
    A
end

###
# FunctionalOperator treats a functional like an operator
###

immutable FunctionalOperator{FT,T} <: BandedOperator{T}
    func::FT
end

FunctionalOperator{T}(func::Functional{T})=FunctionalOperator{typeof(func),T}(func)

Base.convert{T}(::Type{BandedOperator{T}},FO::FunctionalOperator)=FunctionalOperator{typeof(FO.func),T}(FO.func)

bandinds(FO::FunctionalOperator)=0,datalength(FO.func)-1
domainspace(FO::FunctionalOperator)=domainspace(FO.func)
rangespace(FO::FunctionalOperator)=ConstantSpace()

for TYP in (:AnySpace,:UnsetSpace,:FunctionSpace)
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
