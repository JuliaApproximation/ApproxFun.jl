

export MutableOperator



###
# FillMatrix represents the filled-in rows of an almost-bande dmatrix
###

type FillMatrix{T,R}
    bc::Vector{R}         # The boundary rows as functionals
    data::Matrix{T}       # The combination of bcs
    datalength::Int
    numbcs::Int            # The length of bc.  We store this for quicker access, but maybe remove
    padbc::Int             # The bc data needs to padded by this amount to ensure fast backsubstitution,etc.
end




function eye2{BT<:Array}(::Type{BT},n::Integer,m::Integer)
    @assert n==0 || m==0
    Array(BT,n,m)
end
eye2{T}(::Type{T},n::Integer,m::Integer)=eye(T,n,m)
eye2{T}(::Type{T},n::Integer)=eye(T,n,n)

function FillMatrix(bc::Vector,data::Matrix,pf)
    nbc=length(bc)
    data=pad(data,size(data,1)+50,:)  # pad now by 50 to save computational cost later

    ##TODO Maybe better for user to do SavedFunctional?  That way it can be reused
    sfuncs=Array(SavedFunctional{isempty(bc)?Float64:mapreduce(eltype,promote_type,bc)},nbc)

    m=50
    for k=1:nbc
        if isa(bc[k],SavedFunctional)
            sfuncs[k]=bc[k]
        else
            sfuncs[k]=SavedFunctional(bc[k])
        end
        resizedata!(sfuncs[k],m+pf)
    end
    FillMatrix(sfuncs,data,size(data,1),nbc,pf)
end

FillMatrix{T}(::Type{T},bc,pf)=FillMatrix(bc,eye2(T,length(bc)),pf)


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

function unsafe_getindex{T,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    ret = zero(T)

    @simd for m=1:B.numbcs
         @inbounds ret += B.data[k,m]*B.bc[m].data[j]
    end

    ret
end

function getindex{T<:Number,R}(B::FillMatrix{T,R},kr::Range,j::Integer)
    ret = zeros(T,length(kr))

    for m=1:B.numbcs
        @assert j <= B.bc[m].datalength #TODO: temporary for debugging

        bcv = B.bc[m].data[j]
        fd=B.data[kr,m]
        ret += fd*bcv
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


## MutableOperator


type MutableOperator{T,M,R} <: BandedBelowOperator{T}
    op::M                 # The underlying op that is modified
    data::BandedMatrix{T} # Data representing bands
    fill::FillMatrix{T,R}

    shift::Int            # How far down to shift the op

    datalength::Int       # How long data is.  We can't use the array length of data as we double the memory allocation but don't want to fill in

    bandinds::Tuple{Int,Int}   # Encodes the bandrange
end

domainspace(M::MutableOperator)=domainspace(M.op)
rangespace(M::MutableOperator)=rangespace(M.op)



#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableOperator{R<:Functional}(bc::Vector{R},op::BandedOperator)
    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1
    nbc = length(bc)

    br=((bndinds[1]-nbc),(bndindslength-1))
    data = bazeros(op,nbc+100-br[1],:,br)

     # do all columns in the row, +1 for the fill
    fl=FillMatrix(eltype(data),bc,br[end]+1)

    for k=1:nbc,j=columnrange(data,k)
        data[k,j]=fl.bc[k][j]  # initialize data with the boundary rows
    end

    MutableOperator(op,data,fl,nbc,nbc, br)
end

function MutableOperator{T<:Operator}(B::Vector{T})
    bcs = Functional{eltype(eltype(B))}[B[k] for k=1:length(B)-1]

    @assert typeof(B[end]) <: BandedBelowOperator

    MutableOperator(bcs,B[end])
end

MutableOperator{BO<:BandedOperator}(B::BO)=MutableOperator(BO[B])


# for bandrange, we save room for changed entries during Givens
bandinds(B::MutableOperator)=B.bandinds


function Base.getindex{T<:Number,M,R}(B::MutableOperator{T,M,R},kr::UnitRange,jr::UnitRange)
    ret = zeros(T,length(kr),length(jr))


    for k = kr
        if k <= B.datalength
            for j=jr
                ret[k-kr[1]+1,j-jr[1]+1] = B[k,j]
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






function Base.getindex(B::MutableOperator,k::Integer,j::Integer)
    ir = columninds(B,k)::Tuple{Int,Int}
    shift = B.shift

    if k <= B.datalength && j <= ir[end] && ir[1] <= j
        B.data[k,j]
    elseif k <= B.datalength && j > ir[end]
        B.fill[k,j]
    else
        B.op[k-shift,j]##TODO: Slow
    end
end


# getindex!(b::MutableOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
# getindex!(b::MutableOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator,R}(B::MutableOperator{T,M,R},n::Integer)
    resizedata!(B.fill,n)

    if n > B.datalength
        shift=B.shift

        if n > size(B.data,1)
            pad!(B.data,2n,:)
        end

        addentries!(B.op,IndexStride(B.data,shift,0),B.datalength+1-shift:n-shift,:)
        B.datalength = n
    end

    B
end




function Base.setindex!(B::MutableOperator,x,k::Integer,j::Integer)
    resizedata!(B,k)
    #@inbounds
    B.data[k,j] = x
    x
end
