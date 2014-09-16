export FourierDerivativeOperator


## Derivative


function addentries!(D::DerivativeOperator{Complex{Float64},LaurentSpace},A::ShiftArray,kr::Range1)
    d=domain(D)
    m=D.order
    C=2π./(d.b-d.a)*im

    for k=kr
        if isodd(k)
            A[k,0] += (C*(k-1)/2)^m
        else
            A[k,0] += (-C*k/2)^m
        end
    end
    
    A
end


## Multiplication 

addentries!{T}(M::MultiplicationOperator{T,LaurentSpace,LaurentSpace},A,k)=addentries!(LaurentOperator(M.f),A,k)