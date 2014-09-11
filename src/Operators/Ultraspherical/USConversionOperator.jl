export USConversionOperator

## USConversionOperator

type USConversionOperator <: BandedOperator{Float64}
    λ::Int
    domain::IntervalDomain
end

function ConversionOperator{ao,bo}(a::UltrasphericalSpace{ao},b::UltrasphericalSpace{bo})
    @assert domain(a) == domain(b)
    USConversionOperator(ao:bo,domain(a))
end


domainspace(M::USConversionOperator)=UltrasphericalSpace{M.λ-1}(M.domain)
rangespace(M::USConversionOperator)=UltrasphericalSpace{M.λ}(M.domain)

function USConversionOperator(r::Range1,d::IntervalDomain)
    @assert length(r)>1 && r[end] > r[1]

    if length(r)==2
        USConversionOperator(r[2],d)
    else
        USConversionOperator(r[end],d)*USConversionOperator(r[1]:r[end-1],d)
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


function addentries!(C::USConversionOperator,A::ShiftArray,kr::Range1)
    kr = max(kr[1],1):kr[end]

    if C.λ == 1
        one_conversion_addentries!(A,kr)
    else
        conversion_addentries!(C.λ,A,kr)
    end
end


function multiplyentries!(C::USConversionOperator,A::ShiftArray,kr::Range1)
    if C.λ == 1
        one_conversion_multiplyentries!(A,kr)
    else
        conversion_multiplyentries!(C.λ,A,kr)
    end
    
    A
end


bandinds(C::USConversionOperator)=0,2



## spaceconversion

# return the space that has banded ConversionOperator
function conversion_rule{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end


spaceconversion(g::Vector,::UltrasphericalSpace{1},::ChebyshevSpace)=ultraiconversion(g)
spaceconversion(g::Vector,::ChebyshevSpace,::UltrasphericalSpace{1})=ultraconversion(g)
