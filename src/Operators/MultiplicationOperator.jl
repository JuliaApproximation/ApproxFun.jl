export MultiplicationOperator


type ConstantOperator{T<:Number} <: BandedOperator{T}
    c::T
    
    space::Int
end

domainspace(M::ConstantOperator)=M.space
rangespace(M::ConstantOperator)=M.space

bandrange(T::ConstantOperator)=0:0

addentries!(C::ConstantOperator,A::ShiftArray,kr::Range1)=toeplitz_addentries!([.5C.c],A,kr)


type MultiplicationOperator{T<:Number,D<:IntervalDomain} <: BandedOperator{T}
    f::IFun{T,D}
    
    space::Int
end


MultiplicationOperator(c::Number,k::Int)=ConstantOperator(c,k)
MultiplicationOperator(f)=MultiplicationOperator(f,0)



domainspace(M::MultiplicationOperator)=M.space
rangespace(M::MultiplicationOperator)=M.space


function zeromultiplication_addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    toeplitz_addentries!(.5M.f.coefficients,A,kr)
    hankel_addentries!(.5M.f.coefficients,A,max(kr[1],2):kr[end])            
end

function onemultiplication_addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    toeplitz_addentries!(.5M.f.coefficients,A,kr)
    hankel_addentries!(-.5M.f.coefficients[3:end],A,kr)    
end



function usjacobi_addentries!(λ::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,-1]=.5(k-1)/(k-2+λ)
        A[k,1]=.5(k+2λ-1)/(k+λ)
    end
    A
end




usmultiplication_addentries!(λ::Integer,a::IFun,A,kr)=usmultiplication_addentries!(λ,coefficients(a,λ),A,kr)
function usmultiplication_addentries!(λ::Integer,a::Vector,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0]=a[1] 
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1
        J=BandedArray(ShiftArray(zeros(length(jkr),3),2,1-jkr[1]),jkr)
        usjacobi_addentries!(λ,J.data,jkr)
    
        C1=2λ*J;
    
        shiftarray_const_addentries!(C1.data,a[2],A,kr)

        C0=BandedArray(ShiftArray(ones(length(jkr),1),1,1-jkr[1]),jkr)
    
        for k=1:length(a)-2    
            C1,C0=2(k+λ)./(k+1)*J*C1-(k+2λ-1)./(k+1).*C0,C1
            shiftarray_const_addentries!(C1.data,a[k+2],A,kr)    
        end
    end
    
    A
end


function addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    if M.space == 0
        zeromultiplication_addentries!(M,A,kr)
    elseif M.space == 1
        onemultiplication_addentries!(M,A,kr)
    else
        usmultiplication_addentries!(M.space,M.f,A,kr)
    end
end




bandrange(T::MultiplicationOperator)=(1-length(T.f.coefficients):length(T.f.coefficients)-1)
domain(T::MultiplicationOperator)=T.f.domain


