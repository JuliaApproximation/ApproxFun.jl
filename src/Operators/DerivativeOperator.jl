export DerivativeOperator


## DerivativeOperator

type DerivativeOperator{D<:IntervalDomain} <: BandedOperator{Float64}
    order::Range1{Int}
    domain::D
end

DerivativeOperator(k::Integer,d::IntervalDomain)=DerivativeOperator(k-1:k,d)

function addentries!{T<:Interval}(D::DerivativeOperator{T},A::ShiftArray,kr::Range1)
    @assert D.order[1] == 0  ##TODO other orders

    
    μ=D.order[end]
    C=2.^(μ-1).*factorial(μ-1).*(2./(D.domain.b-D.domain.a)).^μ
    
    for k=kr
        A[k,μ] += C*(μ+k-1)
    end
    
    A
end

bandrange(D::DerivativeOperator)=0:(length(D.order)-1)
domainspace(M::DerivativeOperator)=M.order[1]
rangespace(M::DerivativeOperator)=M.order[end]
domain(D::DerivativeOperator)=D.domain


function *(D1::DerivativeOperator,D2::DerivativeOperator)
    @assert D1.domain == D2.domain
    
    DerivativeOperator(D2.order[1]:D2.order[2]+length(D1.order)-1,D1.domain)
end

^(D1::DerivativeOperator,k::Integer)=DerivativeOperator(D1.order[1]:D1.order[1]+k*(length(D1.order)-1),D1.domain)






