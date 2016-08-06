

export MutableOperator



###
# FillMatrix represents the filled-in rows of an almost-bande dmatrix
# data*bc
###

type FillMatrix{T,CO}
    bc::CO                # The boundary rows as a cached operator
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
    if isempty(bc)
        FillMatrix(nothing,data,size(data,1),0,pf)
    elseif length(bc)==1
        FillMatrix(first(bc),data,pf)
    else
        FillMatrix(InterlaceOperator(bc),data,pf)
    end
end

function FillMatrix(bc::Operator,data::Matrix,pf)
    @assert isfinite(size(bc,1))
    nbc=size(data,2)
    data=pad(data,size(data,1)+50,:)  # pad now by 50 to save computational cost later

    sfuncs=cache(bc)
    resizedata!(sfuncs,nbc,size(data,1)+50)
    FillMatrix(sfuncs,data,size(data,1),nbc,pf)
end

FillMatrix{T,OO<:Operator}(::Type{T},bc::Vector{OO},pf) =
    FillMatrix(bc,eye2(T,isempty(bc)?0:mapreduce(op->size(op,1),+,bc)),pf)
FillMatrix{T}(::Type{T},bc::Operator,pf) = FillMatrix(bc,eye2(T,size(bc,1)),pf)



# the bandinds move left as we do row manipulations, so the bandinds[2]
# of F.bc will give the maximum
bandinds(F::FillMatrix) = (-∞,bandinds(F.bc,2))
bandinds{T}(F::FillMatrix{T,Void}) = (0,0)


function getindex{T<:Number,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    ret = zero(T)
    resizedata!(B,k)

    for m=1:B.numbcs
        bcv = B.bc[m,j]
        fd = B.data[k,m]
        ret += fd*bcv
    end

    ret
end

function unsafe_getindex{T,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    ret = zero(T)

    @simd for m=1:B.numbcs
         @inbounds ret += B.data[k,m]*B.bc.data[m,j]
    end

    ret
end

function getindex{T<:Number,R}(B::FillMatrix{T,R},kr::Range,j::Integer)
    ret = zeros(T,length(kr))

    for m=1:B.numbcs
        bcv = B.bc.data[m,j]
        fd=B.data[kr,m]
        ret += fd*bcv
    end

    ret
end

function resizedata!{T}(B::FillMatrix{T},n)
    nbc=B.numbcs
    if nbc>0  && n > B.datalength
        resizedata!(B.bc,:,2n+B.padbc)         ## do all columns in the row, +1 for the fill
        B.data=pad(B.data,2n,:)
        B.datalength=2n
    end
    B
end


## MutableOperator


type MutableOperator{T,M,R} <: Operator{T}
    op::M                 # The underlying op that is modified
    data::BandedMatrix{T} # Data representing bands
    fill::FillMatrix{T,R}

    shift::Int            # How far down to shift the op

    datalength::Int       # How long data is.  We can't use the array length of data as we double the memory allocation but don't want to fill in

    bandinds::Tuple{Int,Int}   # Encodes the bandrange of the banded part
end

domainspace(M::MutableOperator) = domainspace(M.op)

rangespace{T,MM}(M::MutableOperator{T,MM,Void}) = rangespace(M.op)

rangespace(M::MutableOperator) =
    TupleSpace(tuple(spaces(rangespace(M.fill.bc))...,
                     spaces(rangespace(M.op))...))



#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableOperator{R<:Operator}(bc::Vector{R},op::Operator)
    @assert isbanded(op)

    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1
    nbc = isempty(bc)?0:mapreduce(op->size(op,1),+,bc)

    br=((bndinds[1]-nbc),(bndindslength-1))
    data = bzeros(op,nbc+100-br[1],:,br)

     # do all columns in the row, +1 for the fill
    fl=FillMatrix(eltype(data),bc,br[end]+1)

    for k=1:nbc,j=columnrange(data,k)
        data[k,j]=fl.bc[k,j]  # initialize data with the boundary rows
    end

    MutableOperator(op,data,fl,nbc,nbc, br)
end

#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableOperator(bc::Operator,op::Operator)
    @assert isbanded(op)
    nbc = size(bc,1)
    @assert isfinite(nbc)

    bndinds=bandinds(op)
    bndindslength=bndinds[end]-bndinds[1]+1

    br=((bndinds[1]-nbc),(bndindslength-1))
    data = bzeros(op,nbc+100-br[1],:,br)

     # do all columns in the row, +1 for the fill
    fl=FillMatrix(eltype(data),bc,br[end]+1)

    for k=1:nbc,j=columnrange(data,k)
        data[k,j]=fl.bc[k,j]  # initialize data with the boundary rows
    end

    MutableOperator(op,data,fl,nbc,nbc, br)
end


function MutableOperator{T<:Operator}(B::Vector{T})
    bcs = Operator{eltype(eltype(B))}[B[k] for k=1:length(B)-1]

    @assert isinf(size(B[end],1)) && isinf(size(B[end],2))

    MutableOperator(bcs,B[end])
end

MutableOperator{BO<:Operator}(B::BO) = MutableOperator(BO[B])
MutableOperator{T}(B::InterlaceOperator{T,1}) = MutableOperator(B.ops)

function MutableOperator{T}(A::InterlaceOperator{T,2})
    # determine number of bcs and split the operator in two
    nbc=0; while isfinite(size(A.ops[nbc+1],1)) nbc+=1 end
    for k=nbc+2:size(A.ops,2)
        if isfinite(size(A.ops[k],1))
            error("Only the case where all finite operators are at the top is currently supported.")
        end
    end

    #TODO: Remove empty bc case

    bc = nbc==0?Operator{T}[]:InterlaceOperator(A.ops[1:nbc,:])

    MutableOperator(bc,InterlaceOperator(A.ops[nbc+1:end,:]))
end



# for bandrange, we save room for changed entries during Givens
bandinds(B::MutableOperator) = B.bandinds[1],max(B.bandinds[2],bandinds(B.fill,2))


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
    ir = max(k+B.bandinds[1],1),k+B.bandinds[2]
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

function resizedata!{T<:Number,M<:Operator,R}(B::MutableOperator{T,M,R},n::Integer)
    resizedata!(B.fill,n)

    if n > B.datalength
        shift=B.shift

        if n > size(B.data,1)
            pad!(B.data,2n,:)
        end

        kr=B.datalength+1:n
        jr=max(B.datalength+1-B.data.l,1):n+B.data.u
        BLAS.axpy!(1.0,view(B.op,kr-shift,jr),view(B.data,kr,jr))

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
