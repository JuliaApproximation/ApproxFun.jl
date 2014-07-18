export FourierDerivativeOperator



function fourier_derivative_addentries!(m::Integer,d::PeriodicInterval,A::ShiftArray,kr::Range1)
    C=2π./(d.b-d.a)*im

    for k=kr
        A[k,0] += (C*k)^m
    end
    
    A
end




type FourierDerivativeOperator{D<:PeriodicInterval} <: BandedShiftOperator{Complex{Float64}}
    order::Int
    domain::D
end


addentries!(D::FourierDerivativeOperator,A::ShiftArray,kr::Range1)=fourier_derivative_addentries!(D.order,D.domain,A,kr)


bandinds(D::FourierDerivativeOperator)=0,0

domain(D::FourierDerivativeOperator)=D.domain


domainspace(D::FourierDerivativeOperator)=FourierSpace(D.domain)
rangespace(D::FourierDerivativeOperator)=FourierSpace(D.domain)


^(D1::FourierDerivativeOperator,k::Integer)=FourierDerivativeOperator(D1.order*k,D1.domain)