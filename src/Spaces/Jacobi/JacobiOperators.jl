## Evaluation

function Base.getindex(op::Evaluation{Float64,Bool,JacobiSpace},kr::Range1)
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
function Base.getindex(op::Evaluation{Float64,Float64,JacobiSpace},kr::Range1)
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

rangespace(D::Derivative{Float64,JacobiSpace})=JacobiSpace(D.space.a+1,D.space.b+1)



function getdiagonalentry(T::Derivative{Float64,JacobiSpace},k,j)
    if j==0
        0.
    else #j==1
        .5(k+1+T.space.a+T.space.b)
    end
end


## Integral

## Conversion

function Conversion(L::JacobiSpace,M::JacobiSpace)
    @assert M.b>=L.b && M.a>=L.a

    if (M.b == L.b+1 && M.a == L.a) || (M.b == L.b && M.a == L.a+1)
        Conversion{JacobiSpace,JacobiSpace,Float64}(A,B)
    elseif M.b > L.b+1
        Conversion(JacobiSpace(M.a,M.b-1),M)*Conversion(L,JacobiSpace(M.a,M.b-1))    
    elseif M.a > L.a+1
        Conversion(JacobiSpace(M.a-1,M.b),M)*Conversion(L,JacobiSpace(M.a-1,M.b))            
    end
end   

bandinds(C::Conversion{JacobiSpace,JacobiSpace})=(0,1)



function getdiagonalentry(C::Conversion{JacobiSpace,JacobiSpace},k,j)
    T=C.domainspace
    if T.b==C.rangespace.b+1
        if j==0
            k==1?1.:(T.a+T.b+k)/(T.a+T.b+2k-1)
        else
            (T.a+k)./(T.a+T.b+2k+1)
        end    
    elseif T.a==C.rangespace.a+1
        if j==0
            k==1?1.:(T.a+T.b+k)/(T.a+T.b+2k-1)
        else
            -(T.b+k)./(T.a+T.b+2k+1)
        end  
    else
        error("Not implemented")  
    end
end