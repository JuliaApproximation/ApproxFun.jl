export ToeplitzOperator, index



##TODO: Support T<: ShiftVector
type ToeplitzOperator{V<:Vector} <: BandedOperator
    coefficients::V
end

ToeplitzOperator(f::IFun)=ToeplitzOperator(f.coefficients)


index(T::ToeplitzOperator)=length(T.coefficients)

function Base.getindex(T::ToeplitzOperator,kr::Range1,jr::Range1)
    v = T.coefficients

    ret = spzeros(length(kr),length(jr))


    for m=1:length(v)
        for k=(intersect(kr,jr+1-m))
          ret[k-kr[1]+1,k-jr[1]+m] = v[m]         
        end
    
        for k=(intersect(kr+1-m,jr))
          ret[k-kr[1]+m,k-jr[1]+1] += v[m]         
        end    
    end
  
    ret
end




#TODO: assert?
#TODO: we assume A is initialized to 0
#TODO: Allow shift
function copybandedentries(T::ToeplitzOperator,A::ShiftArray,kr::Range1,ksh::Integer,jr::Range1,jsh::Integer)
    v = T.coefficients
    
    
    for k=kr,j=max(jr[1],1-length(v)):min(jr[end],length(v)-1)
        A[k+ksh,j+jsh] = (j ==0) ? 2v[1] : v[abs(j)+1]
    end
    
    A
end


bandrange(T::ToeplitzOperator)=(1-length(T.coefficients):length(T.coefficients)-1)