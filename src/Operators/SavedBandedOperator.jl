

export SavedBandedOperator












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
bandrange(B::SavedBandedOperator)=bandrange(B.op)
datalength(B::SavedBandedOperator)=B.datalength



function addentries!(B::SavedBandedOperator,A::ShiftArray,kr::Range1)       
    resizedata!(B,kr[end])
    br=bandrange(B)

    for k=kr, j=br
        A[k,j]+=B.data[k,j]
    end
    
    A
end

function resizedata!(B::SavedBandedOperator,n::Integer)
    if n > B.datalength
        resize!(B.data,2n,length(bandrange(B)))

        addentries!(B.op,B.data,B.datalength+1:n)
        
        B.datalength = n
    end
    
    B
end


## SavedRowOperator