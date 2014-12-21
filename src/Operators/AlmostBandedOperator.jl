

export MutableAlmostBandedOperator












## MutableAlmostBandedOperator


type MutableAlmostBandedOperator{T,M,R} <: BandedBelowOperator{T}
    bc::Vector{R}         # The boundary rows
    op::M                 # The underlying op that is modified
    data::ShiftArray{T}   # Shifted to encapsolate bandedness  ##TODO: Change to BandedArray
    filldata::Array{T,2}  # The combination of bcs
    
    bcdata::Array{T,2}  # The filled-in data for the boundary rows.  This is a rectangular
    bcfilldata::Array{T,2}
    
    datalength::Int       # How long data is.  We can't use the array length of data as we double the memory allocation but don't want to fill in
    
    bandinds::(Int,Int)   # Encodes the bandrange
    
    numbcs::Int            # The length of bc.  We store this for quicker access, but maybe remove
end

MutableAlmostBandedOperator(bc,op,data,filldata,bcdata,bcfilldata,datalength,bandinds)=MutableAlmostBandedOperator(bc,op,data,filldata,bcdata,bcfilldata,datalength,bandinds,length(bc))

domainspace(M::MutableAlmostBandedOperator)=domainspace(M.op)
rangespace(M::MutableAlmostBandedOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableAlmostBandedOperator{T<:Number,R<:Functional}(bc::Vector{R},op::BandedOperator{T})
    data = ShiftArray(T,index(op))
    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1
    nbc = length(bc)
    
    
    
    bcdata = T[bc[k][j] for k=1:nbc, j=1:bndindslength+nbc-1]   
    bcfilldata = eye(T,nbc)
                
    br=((bndinds[1]-nbc),(bndindslength-1))
    ##TODO complex functionals
    ##TODO Maybe better for user to do SavedFunctional?  That way it can be reused
    sfuncs=Array(SavedFunctional{Float64},length(bc))
    for k=1:length(bc)
        sfuncs[k]=SavedFunctional(bc[k])
    end
    ar0=Array(T,0,nbc)
    MutableAlmostBandedOperator(sfuncs,op,data,ar0,bcdata,bcfilldata,0, br )
end

function MutableAlmostBandedOperator{T<:Operator}(B::Vector{T})
    bcs = Functional[B[k] for k=1:length(B)-1]
    
    @assert typeof(B[end]) <: BandedOperator
    
    MutableAlmostBandedOperator(bcs,B[end])
end

MutableAlmostBandedOperator{BO<:BandedOperator}(B::BO)=MutableAlmostBandedOperator(BO[B])


index(B::MutableAlmostBandedOperator)=index(B.op)::Int


# for bandrange, we save room for changed entries during Givens
bandinds(B::MutableAlmostBandedOperator)=B.bandinds
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
            ir = B.bandinds
            
            for j=max(ir[1]+k,jr[1]):min(ir[end]+k,jr[end])
                ret[k-kr[1]+1,j-jr[1]+1] = B[k,j]
            end
        end
    end
    
    ret
end


##UNSAFE
function fillgetindex{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},k::Integer,j::Integer)
    nbc = B.numbcs
    ret = zero(T)
    
    if k <= nbc
        for m=1:nbc
            @assert j <= B.bc[m].datalength #TODO: temporary for debugging
            @inbounds bcv = B.bc[m].data[j]    
            @inbounds ret += B.bcfilldata[k,m]*bcv
        end
    else
        for m=1:nbc
            @assert j <= B.bc[m].datalength #TODO: temporary for debugging        
            @inbounds bcv = B.bc[m].data[j]
            @inbounds fd=B.filldata[k-nbc,m]
            ret += fd*bcv
        end    
    end
    
    ret
end

function datagetindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    if k <= B.numbcs
        @inbounds ret=B.bcdata[k,j]
    else
        @inbounds ret=B.data[k-B.numbcs,j-k+B.numbcs]
    end
    
    ret
end


function Base.getindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    ir = columninds(B,k)::(Int,Int)
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
#        warn("Should not access op indices")
        B.op[k-nbc,j]##TODO: Slow
    end
end


getindex!(b::MutableAlmostBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::MutableAlmostBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator,R}(B::MutableAlmostBandedOperator{T,M,R},n::Integer)
    if n > B.datalength    
        nbc=B.numbcs
        
        ## TODO: Requires fill with zeros!!  rename resize!
        resize!(B.data,2n,bandrangelength(B))

        if nbc>0      
            for bc in B.bc
                resizedata!(bc,n+B.bandinds[end]+1)         ## do all columns in the row, +1 for the fill
            end
          
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



