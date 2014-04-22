export FourierDerivativeOperator



function fourier_derivative_addentries!(d::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0] += (1.im*k)^d
    end
    
    A
end




type FourierDerivativeOperator <: BandedShiftOperator{Complex{Float64}}
    order::Int
end


addentries!(D::FourierDerivativeOperator,A::ShiftArray,kr::Range1)=fourier_derivative_addentries!(D.order,A,kr)


bandrange(D::FourierDerivativeOperator)=0:0

domain(::FourierDerivativeOperator)=PeriodicInterval()



^(D1::FourierDerivativeOperator,k::Integer)=FourierDerivativeOperator(D1.order*k)