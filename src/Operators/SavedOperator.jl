

export SavedFunctional,SavedBandedOperator





## SavedFunctional

type SavedFunctional{T<:Number,M<:Functional} <: Functional{T}
    op::M
    data::Vector{T}
    datalength::Int
end

function SavedFunctional{T<:Number}(op::Functional{T})
    data = Array(T,0)
    
    SavedFunctional(op,data,0)
end


function Base.getindex(B::SavedFunctional,k::Integer)
    resizedata!(B,k)
    B.data[k] 
end

function Base.getindex(B::SavedFunctional,k::Range)
    resizedata!(B,k[end])
    B.data[k] 
end



function resizedata!(B::SavedFunctional,n::Integer)
    if n > B.datalength
        resize!(B.data,2n)

        B.data[B.datalength+1:n]=B.op[B.datalength+1:n]
        
        B.datalength = n
    end
    
    B
end



## SavedBandedOperator


type SavedBandedOperator{T<:Number,M<:BandedOperator} <: BandedOperator{T}
    op::M
    data::ShiftArray{T}   #Shifted to encapsolate bandedness
    datalength::Int
end






#TODO: index(op) + 1 -> length(bc) + index(op)
function SavedBandedOperator{T<:Number}(op::BandedOperator{T})
    data = ShiftArray(T,index(op))    
    SavedBandedOperator(op,data,0)
end

index(B::SavedBandedOperator)=index(B.op)::Int




domainspace(M::SavedBandedOperator)=domainspace(M.op)
rangespace(M::SavedBandedOperator)=rangespace(M.op)
bandinds(B::SavedBandedOperator)=bandinds(B.op)
datalength(B::SavedBandedOperator)=B.datalength



function addentries!(B::SavedBandedOperator,A::ShiftArray,kr::Range1)       
    resizedata!(B,kr[end])
    br=bandrange(B)

    for k=kr, j=br
        @safastset!(A,@safastget(B.data,k,j)+@safastget(A,k,j) ,k,j)
    end
    
    A
end

function resizedata!(B::SavedBandedOperator,n::Integer)
    if n > B.datalength
        resize!(B.data,2n,bandrangelength(B))

        addentries!(B.op,B.data,B.datalength+1:n)
        
        B.datalength = n
    end
    
    B
end


