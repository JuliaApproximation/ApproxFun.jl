export CompactOperator



type CompactOperator{T<:Number} <: BandedOperator{T}
    matrix::Array{T,2}
end



function hankel_addentries!(v::Vector,A::ShiftArray,kr::Range1)
    for j=1:length(v)
        for k=intersect(kr,1:j)
            A[k,j-2k+1] += v[j]
        end
    end
    
    A
end


addentries!(T::HankelOperator,A::ShiftArray,kr::Range1)=hankel_addentries!(T.coefficients,A,kr)

bandrange(T::HankelOperator)=(1-length(T.coefficients):length(T.coefficients)-1)



## Laurent Operator

type LaurentOperator{T<:Number} <: BandedShiftOperator{T}
    coefficients::ShiftVector{T}
end


addentries!(T::LaurentOperator,A::ShiftArray,kr::Range1)=toeplitz_addentries!(T.coefficients,A,kr)
bandrange(T::LaurentOperator)=(1-length(T.coefficients):length(T.coefficients)-1)

LaurentOperator(f::FFun)=LaurentOperator(f.coefficients)