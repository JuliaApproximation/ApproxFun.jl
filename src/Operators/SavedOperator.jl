

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

domainspace(F::SavedFunctional)=domainspace(F.op)


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
    data::BandedMatrix{T}   #Shifted to encapsolate bandedness
    datalength::Int
    bandinds::(Int,Int)
end






#TODO: index(op) + 1 -> length(bc) + index(op)
function SavedBandedOperator{T<:Number}(op::BandedOperator{T})
    data = BandedMatrix(T,0,:,bandinds(op))   
    SavedBandedOperator(op,data,0,bandinds(op))
end



for OP in (:domain,:domainspace,:rangespace,:(Base.stride))
    @eval $OP(S::SavedBandedOperator)=$OP(S.op)
end

bandinds(B::SavedBandedOperator)=B.bandinds
datalength(B::SavedBandedOperator)=B.datalength



function addentries!(B::SavedBandedOperator,A,kr::Range)       
    resizedata!(B,kr[end])
    
    addentries!(B.data,A,kr)
    
    A
end

function resizedata!(B::SavedBandedOperator,n::Integer)
    if n > B.datalength
        pad!(B.data,2n,:)

        addentries!(B.op,B.data,B.datalength+1:n)
        
        B.datalength = n
    end
    
    B
end




