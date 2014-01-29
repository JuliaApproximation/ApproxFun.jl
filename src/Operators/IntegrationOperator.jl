export IntegrationOperator

## IntegrationOperator

type IntegrationOperator{D<:IntervalDomain} <: BandedOperator
    order::Range1{Int}
    domain::D
end

IntegrationOperator(k::Integer,d::IntervalDomain)=IntegrationOperator(k:k-1,d)

function bandrange(M::IntegrationOperator)
    @assert first(M.order) == 1
    @assert endof(M.order) == 0

    -1:0
end
domainspace(M::IntegrationOperator)=first(M.order)
rangespace(M::IntegrationOperator)=endof(M.order)
domain(Q::IntegrationOperator)=Q.domain

##Very similar to conversion
function addentries!(Q::IntegrationOperator,A::ShiftArray,kr::Range1)
    @assert first(Q.order) == 1
    @assert endof(Q.order) == 0

    for k=kr
        A[k,-1] += .5(Q.domain.b-Q.domain.a)./(k-1)
    end
    
    A    
end


integrate(d::IntervalDomain)=IntegrationOperator(1,d)