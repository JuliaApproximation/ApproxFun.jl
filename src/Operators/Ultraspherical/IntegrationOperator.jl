export IntegrationOperator

## IntegrationOperator

type IntegrationOperator{D<:IntervalDomain} <: BandedOperator{Float64}
    order::Range1{Int}
    domain::D
end

IntegrationOperator(k::Integer,d::IntervalDomain)=IntegrationOperator(k:k-1,d)

function bandinds(M::IntegrationOperator)
    @assert first(M.order) == 1
    @assert last(M.order) == 0

    -1,0
end
domainspace(M::IntegrationOperator)=UltrasphericalSpace(first(M.order),M.domain)
rangespace(M::IntegrationOperator)=UltrasphericalSpace(last(M.order),M.domain)
domain(Q::IntegrationOperator)=Q.domain

##Very similar to conversion
function addentries!(Q::IntegrationOperator,A::ShiftArray,kr::Range1)
    @assert first(Q.order) == 1
    @assert endof(Q.order) == 0

    for k=max(kr[1],2):kr[end]
        A[k,-1] += .5(Q.domain.b-Q.domain.a)./(k-1)
    end
    
    A    
end


