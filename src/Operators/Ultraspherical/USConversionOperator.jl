


function ConversionOperator{a,b}(A::UltrasphericalSpace{a},B::UltrasphericalSpace{b})
    @assert b > a

    if b==a+1
        ConversionOperator{UltrasphericalSpace{a},UltrasphericalSpace{b}}(A,B)
    else
        d=domain(A)
        ConversionOperator(UltrasphericalSpace{b-1}(d),B)*ConversionOperator(A,UltrasphericalSpace{b-1}(d))
    end
end   


function addentries!(M::ConversionOperator{ChebyshevSpace,UltrasphericalSpace{1}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,0] += (k == 1)? 1. : .5
        A[k,2] += -.5        
    end
    
    A    
end

function addentries!{m,λ}(M::ConversionOperator{UltrasphericalSpace{m},UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
    @assert λ==m+1
    for k=kr
        A[k,0] += (λ-1.)./(k - 2. + λ)
        A[k,2] += -(λ-1.)./(k + λ)
    end
    
    A    
end

function multiplyentries!(M::ConversionOperator{ChebyshevSpace,UltrasphericalSpace{1}},A::ShiftArray,kr::Range)
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

function multiplyentries!{m,λ}(M::ConversionOperator{UltrasphericalSpace{m},UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
    @assert λ==m+1
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

bandinds{m,λ}(C::ConversionOperator{UltrasphericalSpace{m},UltrasphericalSpace{λ}})=0,2



## spaceconversion

# return the space that has banded ConversionOperator to the other
function conversion_rule{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end


spaceconversion(g::Vector,::UltrasphericalSpace{1},::ChebyshevSpace)=ultraiconversion(g)
spaceconversion(g::Vector,::ChebyshevSpace,::UltrasphericalSpace{1})=ultraconversion(g)
