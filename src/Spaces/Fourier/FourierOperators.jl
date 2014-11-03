
## Derivative

bandinds(D::Derivative{LaurentSpace})=0,0
bandinds(D::Integral{LaurentSpace})=0,0   


function addentries!(D::Derivative{LaurentSpace},A::ShiftArray,kr::Range1)
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

addentries!(M::Multiplication{LaurentSpace,LaurentSpace},A,k)=addentries!(LaurentOperator(M.f),A,k)

## Converison

function addentries!(C::Conversion{LaurentSpace,FourierSpace},A::ShiftArray,kr::Range)
    for k=kr
        if k==1
            A[k,0]+=1.
        elseif iseven(k)
            A[k,0]+=-1.im
            A[k,1]+=1.im
        else #isodd(k)
            A[k,0]+=1
            A[k,-1]+=1
        end
    end
    A
end
function addentries!(C::Conversion{FourierSpace,LaurentSpace},A::ShiftArray,kr::Range)
    for k=kr
        if k==1
            A[k,0]+=1.
        elseif iseven(k)
            A[k,0]+=0.5im
            A[k,1]+=0.5
        else #isodd(k)
            A[k,0]+=0.5
            A[k,-1]+=-0.5im
        end
    end
    A
end

bandinds(::Conversion{LaurentSpace,FourierSpace})=-1,1
bandinds(::Conversion{FourierSpace,LaurentSpace})=-1,1

function conversion_rule(A::LaurentSpace,B::FourierSpace)
    @assert domainscompatible(A,B)
    B
end



## Real/Imag

for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{TaylorSpace,T}})=FourierSpace(domain(R))
    end
end

bandinds{T}(::RealOperator{ReImSpace{TaylorSpace,T}})=0,2
bandinds{T}(::ImagOperator{ReImSpace{TaylorSpace,T}})=0,1

function addentries!{T}(R::RealOperator{ReImSpace{TaylorSpace,T}},A::ShiftArray,kr::Range)
    for k=kr
        if isodd(k)
            A[k,0]+=1        
        elseif iseven(k)
            A[k,2]+=-1
        end
    end
    A
end

function addentries!{T}(R::ImagOperator{ReImSpace{TaylorSpace,T}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,1]+=1
    end
    A
end





