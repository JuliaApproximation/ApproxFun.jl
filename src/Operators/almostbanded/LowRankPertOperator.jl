struct LowRankPertOperator{OO,LR,T} <: Operator{T}
    op::OO
    pert::LR

    function LowRankPertOperator{OO,LR,T}(a::OO,b::LR) where {OO,LR,T}
        @assert domainspace(a)==domainspace(b)
        @assert rangespace(a)==rangespace(b)
        new{OO,LR,T}(a,b)
    end
end

function LowRankPertOperator(Bin::Operator,Lin::LowRankOperator)
    B,L2=promotedomainspace([Bin,Lin])
    rsp=rangespace(B)  # use rangespace of B because LowRankOperator only
                        # needs convert, and its unlikely that the rangespaces
                        # will be inferred from L
    L=promoterangespace(L2,rsp)

    LowRankPertOperator{typeof(B),typeof(L),promote_type(eltype(L),eltype(B))}(B,L)
end



convert(::Type{Operator{T}},L::LowRankPertOperator) where {T} =
    LowRankPertOperator(Operator{T}(L.op),Operator{T}(L.pert))::Operator{T}
convert(::Type{Operator},V::AbstractVector{OT}) where {OT<:Operator}=LowRankPertOperator(V)


Base.getindex(L::LowRankPertOperator,k::Integer,j::Integer)=L.op[k,j]+L.pert[k,j]


domainspace(L::LowRankPertOperator)=domainspace(L.op)
rangespace(L::LowRankPertOperator)=rangespace(L.op)
datasize(L::LowRankPertOperator,k...)=datasize(L.pert,k...)

israggedbelow(P::LowRankPertOperator) = israggedbelow(P.op) && israggedbelow(P.pert)
bandwidths(P::LowRankPertOperator) = bandwidth(P.op,1)+bandwidth(P.pert,1) , bandwidth(P.op,2)+bandwidth(P.pert,2)

colstop(L::LowRankPertOperator,k::Integer) = max(colstop(L.op,k),colstop(L.pert,k))

for OP in (:promotedomainspace,:promoterangespace)
    @eval $OP(L::LowRankPertOperator,sp::Space)=LowRankPertOperator($OP(L.op,sp),$OP(L.pert,sp))
end



## algebra

for OP in (:+,:-)
    @eval $OP(A::LowRankPertOperator,B::LowRankPertOperator) =
        LowRankPertOperator($OP(A.op,B.op),$OP(A.pert,B.pert))
end

*(L::LowRankPertOperator,f::Fun)=L.op*f+L.pert*f

*(L::LowRankPertOperator,B::LowRankOperator)=L.op*B+L.pert*B
*(B::LowRankOperator,L::LowRankPertOperator)=B*L.op+B*L.pert



*(A::LowRankPertOperator,B::LowRankPertOperator)=A.op*B + A.pert*B

# ambiguity
for TYP in (:TimesOperator,:ZeroOperator,:PlusOperator,:Conversion,:Operator)
    @eval begin
        +(L::LowRankOperator,B::$TYP) = LowRankPertOperator(B,L)
        +(B::$TYP,L::LowRankOperator) = LowRankPertOperator(B,L)

        -(L::LowRankOperator,B::$TYP)=LowRankPertOperator(-B,L)
        -(B::$TYP,L::LowRankOperator)=LowRankPertOperator(B,-L)

        *(L::LowRankPertOperator,B::$TYP)=LowRankPertOperator(L.op*B,L.pert*B)
        *(B::$TYP,L::LowRankPertOperator)=LowRankPertOperator(B*L.op,B*L.pert)
    end
end
