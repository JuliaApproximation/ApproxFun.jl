

export MutableAlmostBandedOperator







## MutableAlmostBandedOperator


type MutableAlmostBandedOperator{T<:Number,M<:BandedOperator,R<:RowOperator} <: BandedBelowOperator{T}
    bc::Vector{R}
    op::M
    data::ShiftArray{T}   #Shifted to encapsolate bandedness
    filldata::Array{T,2}
    
    bcdata::Array{T,2}
    bcfilldata::Array{T,2}
    
    datalength::Int
    
    bandrange::Range1{Int}
end

domainspace(M::MutableAlmostBandedOperator)=domainspace(M.op)
rangespace(M::MutableAlmostBandedOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableAlmostBandedOperator{T<:Number,R<:RowOperator}(bc::Vector{R},op::BandedOperator{T})
    data = ShiftArray(T,index(op))
    
    nbc = length(bc)
    bcdata = T[bc[k][j] for k=1:nbc, j=1:length(bandrange(op))+nbc-1]
    bcfilldata = eye(T,nbc)
                
    br=(bandrange(op)[1]-nbc):(length(bandrange(op))-1)
    MutableAlmostBandedOperator(bc,op,data,Array(T,0,nbc),bcdata,bcfilldata,0, br )
end

function MutableAlmostBandedOperator{T<:Operator}(B::Vector{T})
    bcs = RowOperator[B[k] for k=1:length(B)-1]
    
    @assert typeof(B[end]) <: BandedOperator
    
    MutableAlmostBandedOperator(bcs,B[end])
end


index(B::MutableAlmostBandedOperator)=index(B.op)::Int
numbcs(B::MutableAlmostBandedOperator)=length(B.bc)::Int

# for bandrange, we save room for changed entries during Givens
bandrange(B::MutableAlmostBandedOperator)=B.bandrange
datalength(B::MutableAlmostBandedOperator)=B.datalength




function Base.getindex{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},kr::Range1,jr::Range1)
    ret = spzeros(T,length(kr),length(jr))
    
    
    for k = kr
        if k <= datalength(B) || k <= numbcs(B)
            for j=jr
                ret[k,j] = B[k,j]
            end
        else
            ir = (bandrange(B) + k)::Range1{Int64}
            
            for j=max(ir[1],jr[1]):min(ir[end],jr[end])
                ret[k-kr[1]+1,j-jr[1]+1] = B[k,j]
            end
        end
    end
    
    ret
end

function fillgetindex{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},k::Integer,j::Integer)
    nbc = numbcs(B)
    ret = zero(T)
    
    if k <= nbc
        for m=1:nbc
            bcv = B.bc[m][j]    
            ret += B.bcfilldata[k,m]*bcv
        end
    else
        for m=1:nbc
            bcv = B.bc[m][j]
            fd=B.filldata[k-nbc,m]::T
            ret += fd*bcv
        end    
    end
    
    ret::T
end

function datagetindex{T,M,R}(B::MutableAlmostBandedOperator{T,M,R},k::Integer,j::Integer)  
    nbc = numbcs(B)
    if k <= nbc
        B.bcdata[k,j]::T
    else
        B.data[k-nbc,j-k+nbc]::T
    end
end

function Base.getindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    ir = indexrange(B,k)::Range1{Int64}
    nbc = numbcs(B)
    
    if k <= nbc
        if j <= ir[end]
            B.bcdata[k,j] 
        else
            fillgetindex(B,k,j)
        end
    elseif k-nbc <= datalength(B) && j <= ir[end] && ir[1] <= j
        B.data[k-nbc,j-k+nbc]
    elseif k-nbc <= datalength(B) && j > ir[end]
        fillgetindex(B,k,j)
    else
        B.op[k-nbc,j]##TODO: Slow
    end
end


getindex!(b::MutableAlmostBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::MutableAlmostBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator,R}(B::MutableAlmostBandedOperator{T,M,R},n::Integer)
    l = datalength(B)
    nbc=numbcs(B)::Int
    if n > l
        resize!(B.data,2n,length(B.bandrange))

        if nbc>0        
            newfilldata=zeros(T,2n,nbc)
            newfilldata[1:l,:]=B.filldata[1:l,:]
            B.filldata=newfilldata
        end
        
        addentries!(B.op,B.data,l+1:n)
        
        B.datalength = n
    end
    
    B
end

function Base.setindex!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    nbc = numbcs(B)

    if k<=nbc
        B.bcdata[k,j] = x
    else
        resizedata!(B,k)      
        B.data[k-nbc,j-k+nbc] = x
    end
    x
end

function setfilldata!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    if k<= numbcs(B)
        B.bcfilldata[k,j] = x
    else
        B.filldata[k-numbcs(B),j] = x
    end
end

getfilldata(B::MutableAlmostBandedOperator,k::Integer,j::Integer)=(k<=numbcs(B))?B.bcfilldata[k,j]:B.filldata[k-numbcs(B),j]



