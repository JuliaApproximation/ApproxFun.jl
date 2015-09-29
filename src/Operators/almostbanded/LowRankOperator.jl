export LowRankOperator

immutable LowRankOperator{S<:Space,T} <: AlmostBandedOperator{T}
    U::Vector{Fun{S,T}}
    V::Vector{Functional{T}}

    function LowRankOperator(U::Vector{Fun{S,T}},V::Vector{Functional{T}})
        @assert length(U) == length(V)
        @assert length(U) > 0
        ds=domainspace(first(V))
        for k=2:length(V)
            @assert domainspace(V[k])==ds
        end
        rs=space(first(U))
        for k=2:length(U)
            @assert space(U[k])==rs
        end
        new(U,V)
    end
end



LowRankOperator{S,T}(U::Vector{Fun{S,T}},V::Vector{Functional{T}})=LowRankOperator{S,T}(U,V)
LowRankOperator{S,T1,T2}(U::Vector{Fun{S,T1}},V::Vector{Functional{T2}})=LowRankOperator(convert(Vector{Fun{S,promote_type(T1,T2)}},U),convert(Vector{Fun{M,promote_type(T1,T2)}},V))
LowRankOperator{S,T,FT<:Functional}(U::Vector{Fun{S,T}},V::Vector{FT})=LowRankOperator{S,T}(U,convert(Vector{Functional{eltype(FT)}},V))

function LowRankOperator{FT<:Functional}(B::Vector{FT})
    rsp=TupleSpace(tuple(map(rangespace,B)...,ZeroSpace()))
    LowRankOperator(
        Fun{typeof(rsp),Float64}[Fun([zeros(k-1);1],rsp) for k=1:length(B)],
        B)
end

LowRankOperator(A::Fun,B::Functional)=LowRankOperator([A],[B])



datasize(L::LowRankOperator,k)=k==1?length(L.U):mapreduce(datalength,max,L.V)
datasize(L::LowRankOperator)=datasize(L,1),datasize(L,2)


domainspace(L::LowRankOperator)=domainspace(first(L.V))
rangespace(L::LowRankOperator)=space(first(L.U))
promoterangespace(L::LowRankOperator,sp::Space)=LowRankOperator(map(u->Fun(u,sp),L.U),L.V)
promotedomainspace(L::LowRankOperator,sp::Space)=LowRankOperator(L.U,map(v->promotedomainspace(v,sp),L.V))

function Base.getindex(L::LowRankOperator,k::Integer,j::Integer)
    ret=zero(eltype(L))
    for p=1:length(L.U)
        if k≤length(L.U[p])
            ret+=L.U[p].coefficients[k]*L.V[p][j]
        end
    end
    ret
end



Base.rank(L::LowRankOperator)=length(L.U)


-(L::LowRankOperator)=LowRankOperator(-L.U,L.V)
