export AbstractLowRankOperator, LowRankOperator

abstract AbstractLowRankOperator{T} <: AlmostBandedOperator{T}

immutable LowRankOperator{S<:Space,T} <: AbstractLowRankOperator{T}
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
LowRankOperator{S,T1,T2}(U::Vector{Fun{S,T1}},V::Vector{Functional{T2}})=LowRankOperator(convert(Vector{Fun{S,promote_type(T1,T2)}},U),
                                                                                         convert(Vector{Functional{promote_type(T1,T2)}},V))
LowRankOperator{FF<:Fun,FT<:Functional}(U::Vector{FF},V::Vector{FT})=LowRankOperator(U,convert(Vector{Functional{eltype(FT)}},V))



LowRankOperator(B::AbstractVector,S...)=LowRankOperator(convert(Vector{Functional{Float64}},B),S...)

LowRankOperator(A::Fun,B::Functional)=LowRankOperator([A],[B])



datasize(L::LowRankOperator,k)=k==1?mapreduce(length,max,L.U):mapreduce(datalength,max,L.V)
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

*(L::LowRankOperator,f::Fun)=sum(map((u,v)->u*(v*f),L.U,L.V))

*(L::LowRankOperator,B::BandedOperator)=LowRankOperator(L.U,map(v->v*B,L.V))
*(B::BandedOperator,L::LowRankOperator)=LowRankOperator(map(u->B*u,L.U),L.V)
*(A::LowRankOperator,B::LowRankOperator)=LowRankOperator((A.V*B.U.').'*A.U,B.V)

+(A::LowRankOperator,B::LowRankOperator)=LowRankOperator([A.U;B.U],[A.V;B.V])
-(A::LowRankOperator,B::LowRankOperator)=LowRankOperator([A.U;-B.U],[A.V;B.V])
