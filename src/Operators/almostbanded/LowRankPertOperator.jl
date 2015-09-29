immutable LowRankPertOperator{OO,LR,T} <: AlmostBandedOperator{T}
    op::OO
    pert::LR

    function LowRankPertOperator(a::OO,b::LR)
        @assert domainspace(a)==domainspace(b)
        @assert rangespace(a)==rangespace(b)
        new(a,b)
    end
end

LowRankPertOperator(B::BandedOperator,L::LowRankOperator)=LowRankPertOperator{typeof(B),typeof(L),promote_type(eltype(L),eltype(B))}(B,L)

Base.getindex(L::LowRankPertOperator,k::Integer,j::Integer)=L.op[k,j]+L.pert[k,j]
+(L::LowRankOperator,B::BandedOperator)=LowRankPertOperator(promotespaces([B,L])...)
+(B::BandedOperator,L::LowRankOperator)=LowRankPertOperator(promotespaces([B,L])...)

-(L::LowRankOperator,B::BandedOperator)=LowRankPertOperator(promotespaces([-B,L])...)
-(B::BandedOperator,L::LowRankOperator)=LowRankPertOperator(promotespaces([B,-L])...)

domainspace(L::LowRankPertOperator)=domainspace(L.op)
rangespace(L::LowRankPertOperator)=rangespace(L.op)
datasize(L::LowRankPertOperator,k...)=datasize(L.pert,k...)

for OP in (:promotedomainspace,:promoterangespace)
    @eval $OP(L::LowRankPertOperator,sp::Space)=LowRankPertOperator($OP(L.op,sp),$OP(L.pert,sp))
end

# function MutableOperator(S::LowRankPertOperator)
#     bndinds=bandinds(S.op)
#     dats=datasize(S,1)
#     nbc = 2dats  # need functionals for S.op also
#
#
#     bndindslength=bndinds[end]-bndinds[1]+1
#     br=(bndinds[1],bndindslength+nbc-1)
#
#     data = bazeros(S.op,nbc+100-br[1],br)
#
#      # do all columns in the row, +1 for the fill
#     fl=FillMatrix([S.pert.V;S.op[1:dats,:]],[coefficients(S.pert.U) eye(dats)],br[end]+1)
#
#     for k=1:nbc,j=columnrange(data,k)
#         data[k,j]=S[k,j]  # initialize data with the boundary rows end
#     end
#     MutableOperator(S.op,data,fl,nbc, br)
# end


function MutableOperator{R<:Functional}(bc::Vector{R},S::LowRankPertOperator)
    bndinds=bandinds(S.op)

    dats=datasize(S,1)
    lbc=length(bc)
    shift = lbc+dats


    bndindslength=bndinds[end]-bndinds[1]+1
    br=(bndinds[1]-lbc,bndindslength+dats-1)

    data = bazeros(S.op,shift+100-br[1],br)

    # do all columns in the row, +1 for the fill
    bcdat=eye(shift,lbc)
    lowrdat=[zeros(lbc,dats);coefficients(S.pert.U)]  # add zeros for first bc rows
    opdat=[zeros(lbc,dats);eye(dats)]

    fl=FillMatrix([bc;S.pert.V;S.op[1:dats,:]],[bcdat lowrdat opdat],br[end]+1)


    for k=1:lbc,j=columnrange(data,k)
        data[k,j]=bc[k][j]  # initialize data with the boundary rows
    end

    for k=lbc+1:shift,j=columnrange(data,k)
        data[k,j]=S[k-lbc,j]  # initialize data with the boundary rows end
    end
    M=MutableOperator(S.op[dats+1:end,1:end],data,fl,shift,shift, br)
end
