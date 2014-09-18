## Evaluation

function Base.getindex(op::Evaluation{JacobiSpace,Bool},kr::Range1)
    @assert op.order <= 2
    sp=op.space
    a=sp.a;b=sp.b
    x=op.x
    
    if op.order == 0
        jacobip(kr-1,a,b,x?1.0:-1.0)
    elseif op.order == 1
        @assert !x && b==0 
        Float64[.5*(a+k)*(k-1)*(-1)^k for k=kr]
    elseif op.order == 2
        @assert !x && b==0     
        Float64[-.125*(a+k)*(a+k+1)*(k-2)*(k-1)*(-1)^k for k=kr]
    end
end
function Base.getindex(op::Evaluation{JacobiSpace,Float64},kr::Range1)
    @assert op.order == 0
    jacobip(kr-1,op.space.a,op.space.b,op.x)        
end

## Multiplication

function addentries!{T}(M::Multiplication{T,ChebyshevSpace,JacobiSpace},A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0]=M.f.coefficients[1] 
    end
    
    if length(M.f) > 1
        sp=M.space
        jkr=max(1,kr[1]-length(M.f)+1):kr[end]+length(M.f)-1
        ##TODO: simplify shift array and combine with Ultraspherical
        J=BandedArray(ShiftArray(zeros(length(jkr),3),1-jkr[1],2),jkr)
        addentries!(JacobiRecurrence(sp.a,sp.b).',J.data,jkr)  #Multiplication is transpose
    
        C1=J
    
        shiftarray_const_addentries!(C1.data,M.f.coefficients[2],A,kr)

        C0=BandedArray(ShiftArray(ones(length(jkr),1),1-jkr[1],1),jkr)
    
        for k=1:length(M.f)-2    
            C1,C0=2J*C1-C0,C1
            shiftarray_const_addentries!(C1.data,M.f.coefficients[k+2],A,kr)    
        end
    end
    
    A
end


## Derivative

rangespace(D::Derivative{JacobiSpace})=JacobiSpace(D.space.a+1,D.space.b+1)



function getdiagonalentry(T::Derivative{JacobiSpace},k,j)
    if j==0
        0.
    else #j==1
        .5(k+1+T.space.a+T.space.b)
    end
end


## Integral

## Conversion
# We can only increment by a or b by one, so the following
# multiplies conversion operators to handle otherwise

function Conversion(L::JacobiSpace,M::JacobiSpace)
    @assert (isapprox(M.b,L.b)||M.b>=L.b) && (isapprox(M.a,L.a)||M.a>=L.a)
    
    if (isapprox(M.b,L.b+1) && isapprox(M.a,L.a)) || (isapprox(M.b,L.b) && isapprox(M.a,L.a+1))
        Conversion{JacobiSpace,JacobiSpace,Float64}(L,M)
    elseif M.b > L.b+1
        Conversion(JacobiSpace(M.a,M.b-1),M)*Conversion(L,JacobiSpace(M.a,M.b-1))    
    else  #if M.a >= L.a+1
        Conversion(JacobiSpace(M.a-1,M.b),M)*Conversion(L,JacobiSpace(M.a-1,M.b))            
    end
end   

bandinds(C::Conversion{JacobiSpace,JacobiSpace})=(0,1)



function getdiagonalentry(C::Conversion{JacobiSpace,JacobiSpace},k,j)
    L=C.domainspace
    if L.b+1==C.rangespace.b
        if j==0
            k==1?1.:(L.a+L.b+k)/(L.a+L.b+2k-1)
        else
            (L.a+k)./(L.a+L.b+2k+1)
        end    
    elseif L.a+1==C.rangespace.a
        if j==0
            k==1?1.:(L.a+L.b+k)/(L.a+L.b+2k-1)
        else
            -(L.b+k)./(L.a+L.b+2k+1)
        end  
    else
        error("Not implemented")  
    end
end




# return the space that has banded Conversion to the other
function conversion_rule(A::JacobiSpace,B::JacobiSpace)
    if A.a<=B.a || A.b<=B.b
        A
    else
        B
    end
end



