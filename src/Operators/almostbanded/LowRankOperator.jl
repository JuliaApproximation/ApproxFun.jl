export AbstractLowRankOperator, LowRankOperator

abstract type AbstractLowRankOperator{T} <: Operator{T} end

struct LowRankOperator{S<:Space,T} <: AbstractLowRankOperator{T}
    U::Vector{VFun{S,T}}
    V::Vector{Operator{T}}

    function LowRankOperator{S,T}(U::Vector{VFun{S,T}},V::Vector{Operator{T}}) where {S,T}
        @assert all(isafunctional,V)

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
        new{S,T}(U,V)
    end
end



LowRankOperator(U::Vector{VFun{S,T}},V::Vector{Operator{T}}) where {S,T} = LowRankOperator{S,T}(U,V)
LowRankOperator(U::Vector{VFun{S,T1}},V::Vector{Operator{T2}}) where {S,T1,T2} =
    LowRankOperator(convert(Vector{VFun{S,promote_type(T1,T2)}},U),
                    convert(Vector{Operator{promote_type(T1,T2)}},V))
LowRankOperator(U::Vector{FF},V::Vector{FT}) where {FF<:Fun,FT<:Operator} =
    LowRankOperator(U,convert(Vector{Operator{eltype(FT)}},V))



LowRankOperator(B::AbstractVector,S...) = LowRankOperator(convert(Vector{Operator{Float64}},B),S...)

LowRankOperator(A::Fun,B::Operator) = LowRankOperator([A],[B])


convert(::Type{Operator{T}},L::LowRankOperator{S}) where {S,T} =
    LowRankOperator{S,T}(convert(Vector{VFun{S,T}},L.U),
                         convert(Vector{Operator{T}},L.V))


datasize(L::LowRankOperator,k) =
    k==1?mapreduce(ncoefficients,max,L.U):mapreduce(bandwidth,max,L.V)
datasize(L::LowRankOperator) = datasize(L,1),datasize(L,2)
bandinds(L::LowRankOperator) = 1-datasize(L,1),datasize(L,2)-1

domainspace(L::LowRankOperator) = domainspace(first(L.V))
rangespace(L::LowRankOperator) = space(first(L.U))
promoterangespace(L::LowRankOperator,sp::Space) = LowRankOperator(map(u->Fun(u,sp),L.U),L.V)
promotedomainspace(L::LowRankOperator,sp::Space) = LowRankOperator(L.U,map(v->promotedomainspace(v,sp),L.V))

function Base.getindex(L::LowRankOperator,k::Integer,j::Integer)
    ret=zero(eltype(L))
    for p=1:length(L.U)
        if k≤ncoefficients(L.U[p])
            ret+=L.U[p].coefficients[k]*L.V[p][j]
        end
    end
    ret
end



Base.rank(L::LowRankOperator) = length(L.U)


-(L::LowRankOperator) = LowRankOperator(-L.U,L.V)

*(L::LowRankOperator,f::Fun) = sum(map((u,v)->u*(v*f),L.U,L.V))


*(A::LowRankOperator,B::LowRankOperator) = LowRankOperator((A.V*B.U.').'*A.U,B.V)
# avoid ambiguituy
for TYP in (:TimesOperator,:PlusOperator,:Conversion,:Operator)
    @eval *(L::LowRankOperator,B::$TYP) = LowRankOperator(L.U,map(v->v*B,L.V))
    @eval *(B::$TYP,L::LowRankOperator) = LowRankOperator(map(u->B*u,L.U),L.V)
end

+(A::LowRankOperator,B::LowRankOperator) = LowRankOperator([A.U;B.U],[A.V;B.V])
-(A::LowRankOperator,B::LowRankOperator) = LowRankOperator([A.U;-B.U],[A.V;B.V])
