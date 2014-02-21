export DerivativeOperator



function derivative_addentries!(μ::Integer,d::Interval,A::ShiftArray,kr::Range1)
    C=2.^(μ-1).*factorial(μ-1).*(2./(d.b-d.a)).^μ
    
    for k=kr
        A[k,μ] += C*(μ+k-1)
    end
    
    A
end



## TODO: Unify Defs

type USDerivativeOperator <: BandedOperator{Float64}
    order::Range1{Int}
end


function addentries!(D::USDerivativeOperator,A::ShiftArray,kr::Range1)
    @assert D.order[1] == 0  ##TODO other orders
    derivative_addentries!(D.order[end],Interval(),A,kr)
end

bandrange(D::USDerivativeOperator)=0:(length(D.order)-1)
domainspace(M::USDerivativeOperator)=M.order[1]
rangespace(M::USDerivativeOperator)=M.order[end]


## DerivativeOperator



type DerivativeOperator{D<:Interval} <: BandedOperator{Float64}
    order::Range1{Int}
    domain::D
end

function DerivativeOperator(order::Range1,d::IntervalDomain)
    @assert order[1] == 0 && order[end] <= 2  ##TODO other orders
    
    Mp = Fun(x->tocanonicalD(d,x),d)
    
    if order[end] == 1
        Mp*USDerivativeOperator(0:1)
    elseif order[end] == 2
        (Mp.^2)*USDerivativeOperator(0:2) + diff(Mp)*USDerivativeOperator(0:1)
    end
end

DerivativeOperator(k::Integer,d::IntervalDomain)=DerivativeOperator(k-1:k,d)



## Operator Routine


function addentries!(D::DerivativeOperator,A::ShiftArray,kr::Range1)
    @assert D.order[1] == 0  ##TODO other orders
    derivative_addentries!(D.order[end],D.domain,A,kr)
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






