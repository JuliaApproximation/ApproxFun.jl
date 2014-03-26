export FourierDerivativeOperator



function fourier_derivative_addentries!(d::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0] += (1.im*k)^d
    end
    
    A
end




type FourierDerivativeOperator <: BandedShiftOperator{Float64}
    order::Int
end


addentries!(D::FourierDerivativeOperator,A::ShiftArray,kr::Range1)=fourier_derivative_addentries!(D.order,A,kr)


bandrange(D::FourierDerivativeOperator)=0:0

