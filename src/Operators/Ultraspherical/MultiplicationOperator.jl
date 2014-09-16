function addentries!{T}(M::MultiplicationOperator{T,ChebyshevSpace},A::ShiftArray,kr::Range1)
    cfs=coefficients(M.f)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(.5cfs,A,max(kr[1],2):kr[end])            
end

function addentries!{T}(M::MultiplicationOperator{T,UltrasphericalSpace{1}},A::ShiftArray,kr::Range1)
    cfs=coefficients(M.f)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(-.5cfs[3:end],A,kr)    
end



function usjacobi_addentries!(λ::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,-1]=.5(k-1)/(k-2+λ)
        A[k,1]=.5(k+2λ-1)/(k+λ)
    end
    A
end

function addentries!{T,λ}(M::MultiplicationOperator{T,UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
    a=coefficients(M.f,domainspace(M))
    for k=kr
        A[k,0]=a[1] 
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1
        J=BandedArray(ShiftArray(zeros(length(jkr),3),1-jkr[1],2),jkr)
        usjacobi_addentries!(λ,J.data,jkr)
    
        C1=2λ*J
    
        shiftarray_const_addentries!(C1.data,a[2],A,kr)

        C0=BandedArray(ShiftArray(ones(length(jkr),1),1-jkr[1],1),jkr)
    
        for k=1:length(a)-2    
            C1,C0=2(k+λ)./(k+1)*J*C1-(k+2λ-1)./(k+1).*C0,C1
            shiftarray_const_addentries!(C1.data,a[k+2],A,kr)    
        end
    end
    
    A
end

 

 
