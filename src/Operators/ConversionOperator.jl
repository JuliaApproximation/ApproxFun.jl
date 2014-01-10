export ConversionOperator

## ConversionOperator

type ConversionOperator <: BandedOperator
    λ::Integer
end

function ConversionOperator(r::Range1)
    if length(r)==2
        ConversionOperator(r[1])
    else
        ConversionOperator(r[end]-1)*ConversionOperator(r[1]:r[end-1])
    end
end    

function one_conversion_addentries!(A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0] += (k == 1)? 1. : .5
        A[k,2] += -.5        
    end
    
    A    
end

function conversion_addentries!(λ::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0] += λ./(k - 1 + λ)
        A[k,2] += -λ./(k + λ + 1)
    end
    
    A    
end

function one_conversion_multiplyentries!(A::ShiftArray,kr::Range1)
    cr=columnrange(A)::Range1{Int64}
    
    #We assume here that the extra rows are redundant
    for k=max(2,kr[1]):kr[end]+2,j=cr
        A[k,j] *= .5
    end
    
    #We assume that A has allocated 2 more bandwidth
    for k=kr[1]:kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]::Float64
    end 
end

function conversion_multiplyentries!(λ::Integer,A::ShiftArray,kr::Range1)
    cr=columnrange(A)::Range1{Int64}
    
    λf = 1.λ
    
    #We assume here that the extra rows are redundant
    for k=kr[1]:kr[end]+2,j=cr
        A[k,j] *= λf./(k - 1. + λf)
    end
    
    #We assume that A has allocated 2 more bandwidth
    for k=kr,j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]::Float64
    end 
end


function addentries!(C::ConversionOperator,A::ShiftArray,kr::Range1)
    if C.λ == 0
        one_conversion_addentries!(A,kr)
    else
        conversion_addentries!(C.λ,A,kr)
    end
end


function multiplyentries!(C::ConversionOperator,A::ShiftArray,kr::Range1)
    if C.λ == 0
        one_conversion_multiplyentries!(A,kr)
    else
        conversion_multiplyentries!(C.λ,A,kr)
    end
end


bandrange(C::ConversionOperator)=0:2

