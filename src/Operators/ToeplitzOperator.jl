export ToeplitzOperator, HankelOperator, MultiplicationOperator






type ToeplitzOperator{T<:Number,V<:Union(Vector{T},ShiftVector{T})} <: BandedOperator{T}
    coefficients::V
end

ToeplitzOperator{T<:Number}(V::Vector{T})=ToeplitzOperator{T,typeof(V)}(V)
ToeplitzOperator(f::AbstractFun)=ToeplitzOperator(f.coefficients)




function addentries!{N<:Number,M<:Vector}(T::ToeplitzOperator{N,M},A::ShiftArray,kr::Range1)
    v = T.coefficients
    
    
    for k=kr,j=1-length(v):length(v)-1
        A[k,j] += (j ==0) ? 2v[1] : v[abs(j)+1]
    end
    
    A
end


bandrange{N<:Number,M<:Vector}(T::ToeplitzOperator{N,M})=(1-length(T.coefficients):length(T.coefficients)-1)




function addentries!{N<:Number,M<:ShiftVector}(T::ToeplitzOperator{N,M},A::ShiftArray,kr::Range1)
    v = T.coefficients
    
    
    for k=kr,j=range(v)[1]:range(v)[end]
        A[k,j] += v[j]
    end
    
    A
end


bandrange{N<:Number,M<:ShiftVector}(T::ToeplitzOperator{N,M})=range(T.coefficients)


## Hankel Operator


type HankelOperator{T<:Number} <: BandedOperator{T}
    coefficients::Vector{T}
end

HankelOperator(f::IFun)=HankelOperator(f.coefficients)

function addentries!(T::HankelOperator,A::ShiftArray,kr::Range1)
    v=T.coefficients
  
    for j=1:length(v)
        for k=intersect(kr,1:j)
            A[k,j-2k+1] += v[j]
        end
    end
    
    A
end


bandrange(T::HankelOperator)=(1-length(T.coefficients):length(T.coefficients)-1)


## MultiplicationOperator


type MultiplicationOperator{T<:Number} <: BandedOperator{T}
    T::ToeplitzOperator{T,Vector{T}}
    H::HankelOperator{T}
    
    space::Int
end


MultiplicationOperator(c::Number)=MultiplicationOperator(IFun([1.c]))
MultiplicationOperator(c::Number,k::Int)=MultiplicationOperator(IFun([1.c]),k)
MultiplicationOperator(f::IFun)=MultiplicationOperator(f,0)
function MultiplicationOperator(f::IFun,k::Int)
    if k ==0 
        MultiplicationOperator(ToeplitzOperator(.5f),HankelOperator(.5f),0)
    elseif k ==1
        MultiplicationOperator(ToeplitzOperator(.5f),HankelOperator(-.5f.coefficients[3:end]),1)
    elseif length(f) == 1
        ##TODO: Do a better construction
        MultiplicationOperator(ToeplitzOperator(.5f),HankelOperator(0.f),k)
    else
        error("Higher order multiplication not implemented")
    end
end


domainspace(M::MultiplicationOperator)=M.space
rangespace(M::MultiplicationOperator)=M.space


function zeromultiplication_addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    addentries!(M.T,A,kr)
    addentries!(M.H,A,kr)        
    
    if kr[1] == 1  
        for j=1:length(M.H.coefficients)
            A[1,j-1] -= M.H.coefficients[j]
        end
    end   
    
    A
end

function onemultiplication_addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    addentries!(M.T,A,kr)
    addentries!(M.H,A,kr)        
end

function addentries!(M::MultiplicationOperator,A::ShiftArray,kr::Range1)
    if M.space == 0
        zeromultiplication_addentries!(M,A,kr)
    else
        onemultiplication_addentries!(M,A,kr)
    end
end




bandrange(T::MultiplicationOperator)=(1-length(T.T.coefficients):length(T.T.coefficients)-1)



