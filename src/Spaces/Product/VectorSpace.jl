## VectorSpace{T,S} encodes a space that is a Vector, with coefficients interlaced


immutable VectorDomainSpace{S,T} <: DomainSpace{T}
     space::S
     length::Int
 end

VectorDomainSpace{T}(S::DomainSpace{T},n)=VectorDomainSpace{typeof(S),T}(S,n)

domain(S::VectorDomainSpace)=domain(S.space)
transform(S::VectorDomainSpace,vals::Vector)=transform!(S,hcat(vals...).')


function transform!(S::VectorDomainSpace,M::Array)
    @assert size(M,2)==S.length
    for k=1:size(M,2)
        M[:,k]=transform(S.space,M[:,k])
    end
    vec(M.')
end

Base.vec{S<:DomainSpace,V,T}(f::Fun{VectorDomainSpace{S,V},T})=Fun{S,T}[Fun(f.coefficients[j:f.space.length:end],f.space.space) for j=1:f.space.length]
function devec{S,T}(v::Vector{Fun{S,T}})
    @assert mapreduce(space,isequal,v)
    Fun(vec(coefficients(v).'),VectorDomainSpace(space(first(v)),length(v)))
end

evaluate{V<:VectorDomainSpace,T}(f::Fun{V,T},x)=evaluate(vec(f),x)

for op in (:differentiate,:integrate,:(Base.cumsum))
    @eval $op{V<:VectorDomainSpace}(f::Fun{V})=devec(map($op,vec(f)))
end

