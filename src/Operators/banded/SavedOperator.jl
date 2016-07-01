

export SavedFunctional,SavedBandedOperator





## SavedFunctional

type SavedFunctional{T<:Number,M<:Functional} <: Functional{T}
    op::M
    data::Vector{T}
    datalength::Int
end

SavedFunctional(op::Functional,data)=SavedFunctional(op,data,length(data))
SavedFunctional{T<:Number}(op::Functional{T})=SavedFunctional(op,Array(T,0),0)

for TYP in (:Functional,:Operator)
    @eval Base.convert{T}(::Type{$TYP{T}},S::SavedFunctional)=SavedFunctional(convert(Functional{T},S.op))
end


domainspace(F::SavedFunctional)=domainspace(F.op)
datalength(S::SavedFunctional)=datalength(S.op)

function Base.getindex(B::SavedFunctional,k::Integer)
    resizedata!(B,k)
    B.data[k]
end

function Base.getindex(B::SavedFunctional,k::Range)
    resizedata!(B,k[end]+50)
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
    bandinds::Tuple{Int,Int}
end



# convert needs to throw away calculated data
function Base.convert{BT<:Operator}(::Type{BT},S::SavedBandedOperator)
    T=eltype(BT)
    if isa(S,BT)
        S
    else
        SavedBandedOperator(convert(BandedOperator{T},S.op),
                            convert(BandedMatrix{T},S.data),
                            S.datalength,S.bandinds)
    end
end


#TODO: index(op) + 1 -> length(bc) + index(op)
function SavedBandedOperator{T<:Number}(op::BandedOperator{T})
    data = bzeros(T,0,:,bandinds(op))  # bzeros is needed to allocate top of array
    SavedBandedOperator(op,data,0,bandinds(op))
end



for OP in (:domain,:domainspace,:rangespace,:(Base.stride))
    @eval $OP(S::SavedBandedOperator)=$OP(S.op)
end

bandinds(B::SavedBandedOperator)=B.bandinds
datalength(B::SavedBandedOperator)=B.datalength




function Base.getindex(B::SavedBandedOperator,k::Integer,j::Integer)
    resizedata!(B,k)
    B.data[k,j]
end

function resizedata!(B::SavedBandedOperator,n::Integer)
    if n > B.datalength
        pad!(B.data,2n,:)

        kr=B.datalength+1:n
        jr=max(B.datalength+1-B.data.l,1):n+B.data.u
        BLAS.axpy!(1.0,@compat(view(B.op,kr,jr)),@compat(view(B.data,kr,jr)))

        B.datalength = n
    end

    B
end
