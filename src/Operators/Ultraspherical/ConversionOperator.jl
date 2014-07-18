export ConversionOperator

## ConversionOperator

type ConversionOperator <: BandedOperator{Float64}
    λ::Int
end

domainspace(M::ConversionOperator)=UltrasphericalSpace(M.λ-1)
rangespace(M::ConversionOperator)=UltrasphericalSpace(M.λ)

function ConversionOperator(r::Range1)
    @assert r[end] > r[1]

    if length(r)==2
        ConversionOperator(r[2])
    else
        ConversionOperator(r[end])*ConversionOperator(r[1]:r[end-1])
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
        A[k,0] += (λ-1.)./(k - 2. + λ)
        A[k,2] += -(λ-1.)./(k + λ)
    end
    
    A    
end

function one_conversion_multiplyentries!(A::ShiftArray,kr::Range1)
    cr=columnrange(A)::Range1{Int}
    
    #We assume here that the extra rows are redundant
    for k=max(2,kr[1]):kr[end]+2,j=cr
        A[k,j] *= .5
    end
    
    #We assume that A has allocated 2 more bandwidth
    for k=max(1,kr[1]):kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]
    end 
end

function conversion_multiplyentries!(λ::Integer,A::ShiftArray,kr::Range1)
    cr=columnrange(A)::Range1{Int64}
    
    λf = 1.λ
    
    #We assume here that the extra rows are redundant
    for k=max(kr[1],1):kr[end]+2,j=cr
        A[k,j] *= (λf-1)./(k - 2. + λf)
    end
    
    #We assume that A has allocated 2 more bandwidth
    for k=max(kr[1],1):kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]
    end 
end


function addentries!(C::ConversionOperator,A::ShiftArray,kr::Range1)
    kr = max(kr[1],1):kr[end]

    if C.λ == 1
        one_conversion_addentries!(A,kr)
    else
        conversion_addentries!(C.λ,A,kr)
    end
end


function multiplyentries!(C::ConversionOperator,A::ShiftArray,kr::Range1)
    if C.λ == 1
        one_conversion_multiplyentries!(A,kr)
    else
        conversion_multiplyentries!(C.λ,A,kr)
    end
    
    A
end


bandinds(C::ConversionOperator)=0,2




