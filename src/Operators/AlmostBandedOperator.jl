

export AlmostBandedOperator



###
# FillMatrix represents the filled-in rows of an almost-bande dmatrix
###

type FillMatrix{T,R}
    bc::Vector{R}         # The boundary rows    
    data::Matrix{T}       # The combination of bcs    
    datalength::Int
    numbcs::Int            # The length of bc.  We store this for quicker access, but maybe remove  
    padbc::Int             # The bc data needs to padded by this amount to ensure fast backsubstitution,etc.
end



function FillMatrix{T}(::Type{T},bc,pf)
    nbc=length(bc)

    ##TODO Maybe better for user to do SavedFunctional?  That way it can be reused
    sfuncs=Array(SavedFunctional{isempty(bc)?Float64:mapreduce(eltype,promote_type,bc)},nbc)
    
    m=50
    for k=1:nbc
        sfuncs[k]=SavedFunctional(bc[k])
        resizedata!(sfuncs[k],m+pf)
    end
    ar0=eye(T,m,nbc)  # the first nbc fill in rows are just the bcs
    FillMatrix(sfuncs,ar0,size(ar0,1),nbc,pf)    
end


function getindex{T<:Number,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    ret = zero(T)
    
    for m=1:B.numbcs
        @assert j <= B.bc[m].datalength #TODO: temporary for debugging        

         bcv = B.bc[m].data[j]
        fd=B.data[k,m]
        ret += fd*bcv
    end    
    
    ret
end

function unsafe_getindex{T<:Number,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    ret = zero(T)
    
    @simd for m=1:B.numbcs
         @inbounds ret += B.data[k,m]*B.bc[m].data[j]
    end    
    
    ret
end


function resizedata!{T}(B::FillMatrix{T},n)
    nbc=B.numbcs
    if nbc>0  && n > B.datalength
        for bc in B.bc
            #TODO: Why +10?
            resizedata!(bc,2n+B.padbc)         ## do all columns in the row, +1 for the fill
        end
      
        newfilldata=zeros(T,2n,nbc)
        newfilldata[1:B.datalength,:]=B.data[1:B.datalength,:]
        B.data=newfilldata
        
        B.datalength=2n
    end
    B
end


## AlmostBandedOperator


type AlmostBandedOperator{T,M,R} <: BandedBelowOperator{T}
    op::M                 # The underlying op that is modified
    data::BandedMatrix{T} # Data representing bands
    fill::FillMatrix{T,R}
    
    datalength::Int       # How long data is.  We can't use the array length of data as we double the memory allocation but don't want to fill in
    
    bandinds::(Int,Int)   # Encodes the bandrange
end

domainspace(M::AlmostBandedOperator)=domainspace(M.op)
rangespace(M::AlmostBandedOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function AlmostBandedOperator{T<:Number,R<:Functional}(bc::Vector{R},op::BandedOperator{T})
    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1
    nbc = length(bc)
    
    br=((bndinds[1]-nbc),(bndindslength-1))
    data = bazeros(T,nbc+100-br[1],br)        
    
     # do all columns in the row, +1 for the fill
    fl=FillMatrix(T,bc,br[end]+1) 
    
    for k=1:nbc,j=columnrange(data,k)
        data[k,j]=fl.bc[k][j]  # initialize data with the boundary rows
    end
    
    AlmostBandedOperator(op,data,fl,nbc, br)
end

function AlmostBandedOperator{T<:Operator}(B::Vector{T})
    bcs = Functional[B[k] for k=1:length(B)-1]
    
    @assert typeof(B[end]) <: BandedOperator
    
    AlmostBandedOperator(bcs,B[end])
end

AlmostBandedOperator{BO<:BandedOperator}(B::BO)=AlmostBandedOperator(BO[B])


# for bandrange, we save room for changed entries during Givens
bandinds(B::AlmostBandedOperator)=B.bandinds


function Base.getindex{T<:Number,M,R}(B::AlmostBandedOperator{T,M,R},kr::UnitRange,jr::UnitRange)
    ret = zeros(T,length(kr),length(jr))
    
    
    for k = kr
        if k <= B.datalength
            for j=jr
                ret[k,j] = B[k,j]
            end
        else
            ir = B.bandinds
            
            for j=max(ir[1]+k,jr[1]):min(ir[end]+k,jr[end])
                ret[k-kr[1]+1,j-jr[1]+1] = B[k,j]  #TODO: This is probably slow
            end
        end
    end
    
    ret
end






function Base.getindex(B::AlmostBandedOperator,k::Integer,j::Integer)  
    ir = columninds(B,k)::(Int,Int)
    nbc = B.fill.numbcs
  
    if k <= B.datalength && j <= ir[end] && ir[1] <= j
        B.data[k,j]
    elseif k <= B.datalength && j > ir[end]
        B.fill[k,j]
    else
        B.op[k,j]##TODO: Slow
    end
end


# getindex!(b::AlmostBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
# getindex!(b::AlmostBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator,R}(B::AlmostBandedOperator{T,M,R},n::Integer)
    resizedata!(B.fill,n)

    if n > B.datalength    
        nbc=B.fill.numbcs

        if n > size(B.data,1)
            pad!(B.data,2n)
        end
        
        addentries!(B.op,IndexShift(B.data,-nbc,0),B.datalength+1-nbc:n-nbc)
        B.datalength = n
    end
    
    B
end


function Base.setindex!(B::AlmostBandedOperator,x,k::Integer,j::Integer)
    resizedata!(B,k)      
    #@inbounds 
    B.data[k,j] = x
    x
end


##fast assumes we are inbounds and already resized
# function fastsetindex!(B::AlmostBandedOperator,x,k::Integer,j::Integer)
#     if k<=B.numbcs
#         @inbounds B.bcdata[k,j] = x
#     else
#         ibsetindex!(B.data,x,k-B.numbcs,j-k+B.numbcs)
#     end
#     x
# end



