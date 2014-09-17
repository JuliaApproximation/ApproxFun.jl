
## Derivative


function addentries!(D::Derivative{Complex{Float64},LaurentSpace},A::ShiftArray,kr::Range1)
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

addentries!{T}(M::Multiplication{T,LaurentSpace,LaurentSpace},A,k)=addentries!(LaurentOperator(M.f),A,k)

## Converison

function addentries!(C::Conversion{LaurentSpace,FourierSpace},A::ShiftArray,kr::Range)
    for k=kr
        if k==1
            A[k,0]=1.
        elseif iseven(k)
            A[k,0]=-1.im
            A[k,1]=1.im
        else #isodd(k)
            A[k,0]=1
            A[k,-1]=1
        end
    end
    A
end
function addentries!(C::Conversion{FourierSpace,LaurentSpace},A::ShiftArray,kr::Range)
    for k=kr
        if k==1
            A[k,0]=1.
        elseif iseven(k)
            A[k,0]=0.5im
            A[k,1]=0.5
        else #isodd(k)
            A[k,0]=0.5
            A[k,-1]=-0.5im
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