

export MutableAlmostBandedOperator












## MutableAlmostBandedOperator


type MutableAlmostBandedOperator{T<:Number,M<:BandedOperator,R<:Functional} <: BandedBelowOperator{T}
    bc::Vector{R}
    op::M
    data::ShiftArray{T}   #Shifted to encapsolate bandedness
    filldata::Array{T,2}
    
    bcdata::Array{T,2}
    bcfilldata::Array{T,2}
    
    datalength::Int
    
    bandrange::Range1{Int}
    
    numbcs::Int
end

MutableAlmostBandedOperator(bc,ops...)=MutableAlmostBandedOperator(bc,ops...,length(bc))

domainspace(M::MutableAlmostBandedOperator)=domainspace(M.op)
rangespace(M::MutableAlmostBandedOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableAlmostBandedOperator{T<:Number,R<:Functional}(bc::Vector{R},op::BandedOperator{T})
    data = ShiftArray(T,index(op))
    
    nbc = length(bc)
    bcdata = T[bc[k][j] for k=1:nbc, j=1:length(bandrange(op))+nbc-1]
    bcfilldata = eye(T,nbc)
                
    br=(bandrange(op)[1]-nbc):(length(bandrange(op))-1)
    MutableAlmostBandedOperator(bc,op,data,Array(T,0,nbc),bcdata,bcfilldata,0, br )
end

function MutableAlmostBandedOperator{T<:Operator}(B::Vector{T})
    bcs = Functional[B[k] for k=1:length(B)-1]
    
    @assert typeof(B[end]) <: BandedOperator
    
    MutableAlmostBandedOperator(bcs,B[end])
end


index(B::MutableAlmostBandedOperator)=index(B.op)::Int


# for bandrange, we save room for changed entries during Givens
bandrange(B::MutableAlmostBandedOperator)=B.bandrange
datalength(B::MutableAlmostBandedOperator)=B.datalength



function addentries!(B::MutableAlmostBandedOperator,A::ShiftArray,kr::Range1)
    @assert kr[1] > B.numbcs  #can't write infinite rows 
    
    # We assume that the operator is not filled-in
    
    br=bandrange(B)
    for k=kr[1]:min(kr[end],B.datalength), j=br
        A[k,j]+=B.data[k,j]
    end
    for k=max(B.datalength+1,kr[1]):kr[end], j=br
        A[k,j]+=B.op[k,j]
    end    
    
    A
end



function Base.getindex{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},kr::Range1,jr::Range1)
    ret = spzeros(T,length(kr),length(jr))
    
    
    for k = kr
        if k <= datalength(B) || k <= B.numbcs
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
    nbc = B.numbcs
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

function datagetindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    if k <= B.numbcs
        B.bcdata[k,j]
    else
        B.data[k-B.numbcs,j-k+B.numbcs]
    end
end


function Base.getindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    ir = indexrange(B,k)::Range1{Int64}
    nbc = B.numbcs
    
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
    if n > B.datalength
        nbc=B.numbcs
        resize!(B.data,2n,length(B.bandrange))

        if nbc>0        
            newfilldata=zeros(T,2n,nbc)
            newfilldata[1:B.datalength,:]=B.filldata[1:B.datalength,:]
            B.filldata=newfilldata
        end
        
        addentries!(B.op,B.data,B.datalength+1:n)
        
        B.datalength = n
    end
    
    B
end


function Base.setindex!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    if k<=B.numbcs
        B.bcdata[k,j] = x
    else
        resizedata!(B,k)      
        B.data[k-B.numbcs,j-k+B.numbcs] = x
    end
    x
end


##fast assumes we are inbounds and already resized
function fastsetindex!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    if k<=B.numbcs
        @inbounds B.bcdata[k,j] = x
    else
        @inbounds B.data[k-B.numbcs,j-k+B.numbcs] = x
    end
    x
end

function setfilldata!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    if k<= B.numbcs
        B.bcfilldata[k,j] = x
    else
        B.filldata[k-B.numbcs,j] = x
    end
end

getfilldata(B::MutableAlmostBandedOperator,k::Integer,j::Integer)=(k<=B.numbcs)?B.bcfilldata[k,j]:B.filldata[k-B.numbcs,j]



