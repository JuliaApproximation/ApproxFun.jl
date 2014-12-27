export ToeplitzOperator, HankelOperator, LaurentOperator






type ToeplitzOperator{T<:Number,V<:Union(Vector,ShiftVector)} <: BandedOperator{T}
    coefficients::V
end

ToeplitzOperator{T<:Number}(V::Vector{T})=ToeplitzOperator{T,typeof(V)}(V)
ToeplitzOperator{T<:Number}(V::ShiftVector{T})=ToeplitzOperator{T,typeof(V)}(V)
ToeplitzOperator{T,D}(f::Fun{D,T})=ToeplitzOperator(f.coefficients)




function toeplitz_addentries!(v::Vector,A,kr::UnitRange)    
    if !isempty(v)
        v1=v[1]
        for k=kr
            A[k,k]+=v1
        end
    
        for j=2:length(v)
            vj=v[j]
            for k = kr
                A[k,k+j-1]+=vj
            end
            for k = max(kr[1],j):kr[end]
                A[k,k-j+1]+=vj
            end            
        end    
    end
    A
end

function toeplitz_addentries!(v::Vector,A::BandedMatrix,kr::UnitRange)    
    if !isempty(v)
        @inbounds v1=v[1]
        @simd for k=kr
            @inbounds A.data[A.l+1,k]+=v1
        end
    
        for j=2:length(v)
            @inbounds vj=v[j]
            @simd for k = kr
                @inbounds A.data[j+A.l,k] +=vj
            end
            @simd for k = max(kr[1],j):kr[end]
                @inbounds A.data[2-j+A.l,k]+=vj
            end            
        end    
    end
    A
end

function toeplitz_addentries!(v::ShiftVector,A,kr::Range)    
    for k=kr,j=max(range(v)[1],1-k):range(v)[end]
        A[k,k+j] += v[j]
    end
    
    A
end



addentries!(T::ToeplitzOperator,A,kr::Range)=toeplitz_addentries!(T.coefficients,A,kr)



bandinds{N<:Number,M<:Vector}(T::ToeplitzOperator{N,M})=(1-length(T.coefficients),length(T.coefficients)-1)
bandinds{N<:Number,M<:ShiftVector}(T::ToeplitzOperator{N,M})=firstindex(T.coefficients),lastindex(T.coefficients)


## Hankel Operator


type HankelOperator{T<:Number} <: BandedOperator{T}
    coefficients::Vector{T}
end

HankelOperator(f::Fun)=HankelOperator(f.coefficients)

function hankel_addentries!(v::Vector,A,kr::Range)
    for j=1:length(v)
        for k=max(first(kr),1):min(last(kr),j)
            if j + 1 >= k+1
                A[k,j-k+1] += v[j]
            end
        end
    end
    
    A
end


addentries!(T::HankelOperator,A,kr::Range1)=hankel_addentries!(T.coefficients,A,kr)

bandinds(T::HankelOperator)=(1-length(T.coefficients),length(T.coefficients)-1)



## Laurent Operator

type LaurentOperator{T<:Number} <: BandedOperator{T}
    coefficients::ShiftVector{T}
end



function shiftrowrange(kr::Range)
    if isodd(kr[1]) && isodd(kr[end])
        poskr=div(kr[1]-1,2):div(kr[end]-1,2)
        negkr=-div(kr[end]-1,2):-div(kr[1]+1,2)
    elseif isodd(kr[1]) # && iseven(kr[end])
        poskr=div(kr[1]-1,2):div(kr[end],2)-1
        negkr=-div(kr[end],2):-div(kr[1]+1,2)        
    elseif isodd(kr[end]) # && iseven(kr[1])
        poskr=div(kr[1],2):div(kr[end]-1,2)
        negkr=-div(kr[end]-1,2):-div(kr[1],2)    
    else # iseven(kr[end]) && iseven(kr[1])
        poskr=div(kr[1],2):div(kr[end],2)-1
        negkr=-div(kr[end],2):-div(kr[1],2)        
    end
    negkr,poskr
end

# function laurent_addentries!(v::Vector,A,kr::Range)    
#     negkr,poskr=shiftrowrange(kr)
#     br=1-length(v):length(v)-1
# 
#     for k=poskr,j=br
#         if j≥0 # && k≥0
#             A[2k+1,2j+1] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         else # && k≥0
#             A[2k+1,-2j] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         end
#     end
#     
#     for k=negkr,j=br
#         if j≥0 # && k <0
#             A[-2k,2j+1] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         else # j<0 && k<0
#             A[-2k,-2j] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         end
#     end    
#     
#     A
# end

function laurent_addentries!(v::ShiftVector,A,kr::Range)    
    negkr,poskr=shiftrowrange(kr)
    br=length(v)==0?(0:0):(firstindex(v):lastindex(v))

    for k=poskr,j=br
        if k+j≥0 # && k≥0
            A[2k+1,2k+1+2j] += v[j]
        else # k+j<0 && k≥0
            A[2k+1,2k-2j] += v[j]            
        end
    end
    
    for k=negkr,j=br
        if k+j<0 # && k <0
            A[-2k,-2k-2j] += v[j]
        else # +kj<0 && k<0
            A[-2k,1-2k-2j] += v[j] 
        end
    end    
    
    A
end




addentries!(T::LaurentOperator,A,kr::Range)=laurent_addentries!(T.coefficients,A,kr)

shiftbandinds(T::LaurentOperator)=length(T.coefficients)==0?(0,0):(firstindex(T.coefficients),lastindex(T.coefficients))
function bandinds(T::LaurentOperator)
    sbi=shiftbandinds(T)
    min(2sbi[1],-2sbi[end]),max(2sbi[end],-2sbi[1])
end

