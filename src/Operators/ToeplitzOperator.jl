export ToeplitzOperator, HankelOperator





type ToeplitzOperator{V<:Union(Vector,ShiftVector)} <: BandedOperator
    coefficients::V
end

ToeplitzOperator(f::AbstractFun)=ToeplitzOperator(f.coefficients)
index{M<:Vector}(T::ToeplitzOperator{M})=length(T.coefficients)
index{M<:ShiftVector}(T::ToeplitzOperator{M})=T.coefficients.index



#TODO: assert?
#TODO: we assume A is initialized to 0
function addentries!{M<:Vector}(T::ToeplitzOperator{M},A::ShiftArray,kr::Range1,ksh::Integer,jsh::Integer)
    v = T.coefficients
    
    
    for k=kr,j=1-length(v):length(v)-1
        A[k+ksh,j+jsh] = (j ==0) ? 2v[1] : v[abs(j)+1]
    end
    
    A
end


bandrange{M<:Vector}(T::ToeplitzOperator{M})=(1-length(T.coefficients):length(T.coefficients)-1)



function addentries!{M<:ShiftVector}(T::ToeplitzOperator{M},A::ShiftArray,kr::Range1,ksh::Integer,jsh::Integer)
    v = T.coefficients
    
    
    for k=kr,j=range(v)[1]:range(v)[end]
        A[k+ksh,j+jsh] = v[j]
    end
    
    A
end


bandrange{M<:ShiftVector}(T::ToeplitzOperator{M})=range(T.coefficients)


## Hankel Operator



HankelOperator(f::IFun)=HankelOperator(f.coefficients)
index(T::HankelOperator)=length(T.coefficients)




#TODO: assert?
#TODO: we assume A is initialized to 0
function addentries!(T::HankelOperator,A::ShiftArray,kr::Range1,ksh::Integer,jsh::Integer)
    v=T.coefficients
  
    for j=1:length(v)
        for k=intersect(kr,1:j)
      A[k+ksh,j-2k+1+jsh]=v[j]
        end
    end
end


bandrange(T::HankelOperator)=(1-length(T.coefficients):length(T.coefficients)-1)




# type HankelOperator <: BandedOperator
#     coefficients::Vector
# end
# 
# HankelOperator(f::IFun)=HankelOperator(f.coefficients)
# index(T::HankelOperator)=length(T.coefficients)
# 
# 
# 
# 
# #TODO: assert?
# #TODO: we assume A is initialized to 0
# function addentries!(T::HankelOperator,A::ShiftArray,kr::Range1,ksh::Integer,jr::Range1,jsh::Integer)
#     v = T.coefficients
#     
#     
#     for j=intersect(1:length(v),kr,jr)
#         for k=j-n+1,1):min(j,n)
#             ret[k+ksh,j+jsh-k+1] = v[j]        
#         end
#     end    
#     
#     
#     for k=kr,j=max(jr[1],1-length(v)):min(jr[end],length(v)-1)
#         A[k+ksh,j+jsh] = (j ==0) ? 2v[1] : v[abs(j)+1]
#     end
#     
#     A
# end
# 
# 
# bandrange{M<:Vector}(T::HankelOperator{M})=(1-length(T.coefficients):length(T.coefficients)-1)
