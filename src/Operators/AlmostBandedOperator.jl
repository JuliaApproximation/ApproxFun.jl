

export AlmostBandedOperator












## AlmostBandedOperator


type AlmostBandedOperator{T,M,R} <: BandedBelowOperator{T}
    bc::Vector{R}         # The boundary rows
    op::M                 # The underlying op that is modified
    data::ShiftMatrix{T}  # Data representing bands
    filldata::Matrix{T}   # The combination of bcs
    
    bcdata::Matrix{T}     # The filled-in data for the boundary rows.  This is a rectangular
    bcfilldata::Matrix{T}
    
    datalength::Int       # How long data is.  We can't use the array length of data as we double the memory allocation but don't want to fill in
    
    bandinds::(Int,Int)   # Encodes the bandrange
    
    numbcs::Int            # The length of bc.  We store this for quicker access, but maybe remove
end

AlmostBandedOperator(bc,op,data,filldata,bcdata,bcfilldata,datalength,bandinds)=AlmostBandedOperator(bc,op,data,filldata,bcdata,bcfilldata,datalength,bandinds,length(bc))

domainspace(M::AlmostBandedOperator)=domainspace(M.op)
rangespace(M::AlmostBandedOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function AlmostBandedOperator{T<:Number,R<:Functional}(bc::Vector{R},op::BandedOperator{T})
    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1
    nbc = length(bc)
    
    
    
    bcdata = T[bc[k][j] for k=1:nbc, j=1:bndindslength+nbc-1]   
    bcfilldata = eye(T,nbc)
                
    br=((bndinds[1]-nbc),(bndindslength-1))
    ##TODO Maybe better for user to do SavedFunctional?  That way it can be reused
    sfuncs=Array(SavedFunctional{isempty(bc)?Float64:mapreduce(eltype,promote_type,bc)},length(bc))
    for k=1:length(bc)
        sfuncs[k]=SavedFunctional(bc[k])
    end
    ar0=Array(T,0,nbc)
    
    data = ShiftMatrix(T,0,(br[1]+nbc,br[end]+nbc))    
    AlmostBandedOperator(sfuncs,op,data,ar0,bcdata,bcfilldata,0, br )
end

function AlmostBandedOperator{T<:Operator}(B::Vector{T})
    bcs = Functional[B[k] for k=1:length(B)-1]
    
    @assert typeof(B[end]) <: BandedOperator
    
    AlmostBandedOperator(bcs,B[end])
end

AlmostBandedOperator{BO<:BandedOperator}(B::BO)=AlmostBandedOperator(BO[B])


index(B::AlmostBandedOperator)=index(B.op)::Int


# for bandrange, we save room for changed entries during Givens
bandinds(B::AlmostBandedOperator)=B.bandinds
datalength(B::AlmostBandedOperator)=B.datalength


function addentries!(B::AlmostBandedOperator,A,kr::Range1)
    @assert kr[1] > B.numbcs  #can't write infinite rows 
    
    # We assume that the operator is not filled-in
    
    br=bandrange(B)
    for k=kr[1]:min(kr[end],B.datalength), j=br
        A[k,j]+=B.data[k,j]
    end
    
    addentries!(B.op,A,max(B.datalength+1,kr[1]):kr[end])
end



function Base.getindex{T<:Number,M,R}(B::AlmostBandedOperator{T,M,R},kr::Range1,jr::Range1)
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
function fillgetindex{T<:Number,M,R}(B::AlmostBandedOperator{T,M,R},k::Integer,j::Integer)
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

function datagetindex(B::AlmostBandedOperator,k::Integer,j::Integer)  
    if k <= B.numbcs
        @inbounds ret=B.bcdata[k,j]
    else
        @inbounds ret=B.data[k-B.numbcs,j-k+B.numbcs]
    end
    
    ret
end


function Base.getindex(B::AlmostBandedOperator,k::Integer,j::Integer)  
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


getindex!(b::AlmostBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::AlmostBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator,R}(B::AlmostBandedOperator{T,M,R},n::Integer)
    if n > B.datalength    
        nbc=B.numbcs
        if n > size(B.data,1)
            pad!(B.data,2n)
        end

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


function Base.setindex!(B::AlmostBandedOperator,x,k::Integer,j::Integer)
    if k<=B.numbcs
        B.bcdata[k,j] = x
    else
        resizedata!(B,k)      
        @inbounds B.data[k-B.numbcs,j-k+B.numbcs] = x
    end
    x
end


##fast assumes we are inbounds and already resized
function fastsetindex!(B::AlmostBandedOperator,x,k::Integer,j::Integer)
    if k<=B.numbcs
        @inbounds B.bcdata[k,j] = x
    else
        ibsetindex!(B.data,x,k-B.numbcs,j-k+B.numbcs)
    end
    x
end

function setfilldata!(B::AlmostBandedOperator,x,k::Integer,j::Integer)
    if k<= B.numbcs
        B.bcfilldata[k,j] = x
    else
        B.filldata[k-B.numbcs,j] = x
    end
end

getfilldata(B::AlmostBandedOperator,k::Integer,j::Integer)=(k<=B.numbcs)?B.bcfilldata[k,j]:B.filldata[k-B.numbcs,j]



