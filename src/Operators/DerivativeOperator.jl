export DerivativeOperator


## DerivativeOperator

type DerivativeOperator{D<:IntervalDomain} <: BandedOperator
    order::Range1{Int}
    domain::D
end

DerivativeOperator(k::Integer,d::IntervalDomain)=DerivativeOperator(k-1:k,d)

function addentries!(D::DerivativeOperator,A::ShiftArray,kr::Range1)
    @assert D.order[1] == 0  ##TODO other orders

    
    μ=D.order[end]
    C=2.^(μ-1).*factorial(μ-1).*(2./(D.domain.b-D.domain.a)).^μ
    
    for k=kr
        A[k,μ] += C*(μ+k-1)
    end
    
    A
end

bandrange(D::DerivativeOperator)=0:length(D.order)
domainspace(M::DerivativeOperator)=M.order[1]
rangespace(M::DerivativeOperator)=M.order[end]




Base.diff(d::IntervalDomain,μ::Integer)=DerivativeOperator(0:μ,d)
Base.diff(d::IntervalDomain)=Base.diff(d,1)
Base.eye(d::IntervalDomain)=MultiplicationOperator(IFun([1.],d))




