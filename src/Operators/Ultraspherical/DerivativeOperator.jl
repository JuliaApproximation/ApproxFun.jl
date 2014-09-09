export DerivativeOperator



function derivative_addentries!(order::Range1,d::Interval,A::ShiftArray,kr::Range1)
    m=length(order)-1


    if order[1] == 0
        C=2.^(m-1).*factorial(order[end]-1)*(2./(d.b-d.a)).^m    
        for k=kr
            A[k,m] += C*(order[end]+k-1)
        end
    else
        C=2.^m.*factorial(order[end]-1)./factorial(order[1]-1)*(2./(d.b-d.a)).^m        
        for k=kr        
            A[k,m] += C
        end
    end
    
    A
end



## TODO: Unify Defs

type USDerivativeOperator <: BandedOperator{Float64}
    order::Range1{Int}
end


addentries!(D::USDerivativeOperator,A::ShiftArray,kr::Range1)=derivative_addentries!(D.order,Interval(),A,max(kr[1],1):kr[end])


bandinds(D::USDerivativeOperator)=0,(length(D.order)-1)
domainspace(M::USDerivativeOperator)=UltrasphericalSpace{M.order[1]}()
rangespace(M::USDerivativeOperator)=UltrasphericalSpace{M.order[end]}()


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
DerivativeOperator(d::IntervalDomain)=DerivativeOperator(1,d)


## Operator Routine


addentries!(D::DerivativeOperator,A::ShiftArray,kr::Range1)=derivative_addentries!(D.order,D.domain,A,  max(kr[1],1):kr[end])





bandinds(D::DerivativeOperator)=0,(length(D.order)-1)
domainspace(M::DerivativeOperator)=UltrasphericalSpace{M.order[1]}(M.domain)
rangespace(M::DerivativeOperator)=UltrasphericalSpace{M.order[end]}(M.domain)
domain(D::DerivativeOperator)=D.domain



#promoting domain space is allowed to change range space
promotedomainspace{sporder}(D::DerivativeOperator,sp::UltrasphericalSpace{sporder})=DerivativeOperator(D.order -D.order[1] + sporder,D.domain)


## simplify higher order derivatives
function *(D1::DerivativeOperator,D2::DerivativeOperator)
    @assert D1.domain == D2.domain
    
    DerivativeOperator(D2.order[1]:D2.order[end]+length(D1.order)-1,D1.domain)
end

^(D1::DerivativeOperator,k::Integer)=DerivativeOperator(D1.order[1]:D1.order[1]+k*(length(D1.order)-1),D1.domain)



