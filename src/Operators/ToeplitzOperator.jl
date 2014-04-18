export ToeplitzOperator, HankelOperator, LaurentOperator






type ToeplitzOperator{T<:Number,V<:Union(Vector{T},ShiftVector{T})} <: BandedOperator{T}
    coefficients::V
end

ToeplitzOperator{T<:Number}(V::Vector{T})=ToeplitzOperator{T,typeof(V)}(V)
ToeplitzOperator{T<:Number}(V::ShiftVector{T})=ToeplitzOperator{T,typeof(V)}(V)
ToeplitzOperator(f::AbstractFun)=ToeplitzOperator(f.coefficients)


function toeplitz_addentries!(v::Vector,A::ShiftArray,kr::Range1)    
    for k=kr,j=1-length(v):length(v)-1
        A[k,j] += (j ==0) ? 2v[1] : v[abs(j)+1]
    end
    
    A
end

function toeplitz_addentries!(v::ShiftVector,A::ShiftArray,kr::Range1)    
    for k=kr,j=range(v)[1]:range(v)[end]
        A[k,j] += v[j]
    end
    
    A
end


addentries!(T::ToeplitzOperator,A::ShiftArray,kr::Range1)=toeplitz_addentries!(T.coefficients,A,kr)



bandrange{N<:Number,M<:Vector}(T::ToeplitzOperator{N,M})=(1-length(T.coefficients):length(T.coefficients)-1)
bandrange{N<:Number,M<:ShiftVector}(T::ToeplitzOperator{N,M})=range(T.coefficients)


## Hankel Operator


type HankelOperator{T<:Number} <: BandedOperator{T}
    coefficients::Vector{T}
end

HankelOperator(f::IFun)=HankelOperator(f.coefficients)

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
bandrange(T::LaurentOperator)=firstindex(T.coefficients):lastindex(T.coefficients)

LaurentOperator(f::FFun)=LaurentOperator(f.coefficients)