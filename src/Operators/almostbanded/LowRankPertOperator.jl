immutable LowRankPertOperator{OO,LR,T} <: AlmostBandedOperator{T}
    op::OO
    pert::LR

    function LowRankPertOperator(a::OO,b::LR)
        @assert domainspace(a)==domainspace(b)
        @assert rangespace(a)==rangespace(b)
        new(a,b)
    end
end

function LowRankPertOperator(Bin::BandedOperator,Lin::LowRankOperator)
    B,L=promotespaces(Bin,Lin)
    LowRankPertOperator{typeof(B),typeof(L),promote_type(eltype(L),eltype(B))}(B,L)
end

# convert a vector of functionals and an operator to a LowRnakPertOperator
function LowRankPertOperator{OT<:Operator}(A::Vector{OT})
    A=promotedomainspace(A)
    for k=1:length(A)-1
        @assert isa(A[k],Functional)
    end
    @assert isa(A[end],BandedOperator)
    L=LowRankOperator(A[1:end-1])
    BB=A[end]
    S=SpaceOperator(StrideOperator(BB,length(A)-1,0),domainspace(BB),
                        TupleSpace(map(rangespace,A)))
    L+S
end

for TYP in (:BandedBelowOperator,:AlmostBandedOperator,:Operator)
    @eval Base.convert{OT<:Operator}(::Type{$TYP},V::Vector{OT})=LowRankPertOperator(V)
end


Base.getindex(L::LowRankPertOperator,k::Integer,j::Integer)=L.op[k,j]+L.pert[k,j]


domainspace(L::LowRankPertOperator)=domainspace(L.op)
rangespace(L::LowRankPertOperator)=rangespace(L.op)
datasize(L::LowRankPertOperator,k...)=datasize(L.pert,k...)

for OP in (:promotedomainspace,:promoterangespace)
    @eval $OP(L::LowRankPertOperator,sp::Space)=LowRankPertOperator($OP(L.op,sp),$OP(L.pert,sp))
end








function MutableOperator{R<:Functional}(bc::Vector{R},S::LowRankPertOperator)
    bndinds=bandinds(S.op)

    dats= datasize(S,1)
    lbc=length(bc)
    shift = lbc+dats
    r=rank(S.pert)

    bndindslength=bndinds[end]-bndinds[1]+1
    br=(bndinds[1]-lbc,bndindslength+dats-1)

    data = bazeros(S.op,shift+100-br[1],:,br)

    # do all columns in the row, +1 for the fill
    bcdat=eye(shift,lbc)
    lowrdat=[zeros(lbc,r);coefficients(S.pert.U)]  # add zeros for first bc rows
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






## algebra

+(L::LowRankOperator,B::BandedOperator)=LowRankPertOperator(B,L)
+(B::BandedOperator,L::LowRankOperator)=LowRankPertOperator(B,L)

-(L::LowRankOperator,B::BandedOperator)=LowRankPertOperator(-B,L)
-(B::BandedOperator,L::LowRankOperator)=LowRankPertOperator(B,-L)

for OP in (:+,:-)
    @eval $OP(A::LowRankPertOperator,B::LowRankPertOperator)=LowRankPertOperator($OP(A.op,B.op),$OP(A.pert,B.pert))
end

*(L::LowRankPertOperator,f::Fun)=L.op*f+L.pert*f
*(L::LowRankPertOperator,B::BandedOperator)=LowRankPertOperator(L.op*B,L.pert*B)
*(B::BandedOperator,L::LowRankPertOperator)=LowRankPertOperator(B*L.op,B*L.pert)

*(L::LowRankPertOperator,B::LowRankOperator)=L.op*B+L.pert*B
*(B::LowRankOperator,L::LowRankPertOperator)=B*L.op+B*L.pert



*(A::LowRankPertOperator,B::LowRankPertOperator)=A.op*B + A.pert*B
