
export ProductDomain


immutable ProductDomain{D,T,dim} <: Domain{T,dim}
    domains::D
end

ProductDomain(d::Tuple)=ProductDomain{typeof(d),mapreduce(eltype,promote_type,d),mapreduce(ndims,+,d)}(d)

fromcanonical(d::BivariateDomain,x::Tuple)=fromcanonical(d,x...)
tocanonical(d::BivariateDomain,x::Tuple)=tocanonical(d,x...)

# product domains are their own canonical domain
for OP in (:fromcanonical,:tocanonical)
    @eval $OP(::ProductDomain,x,y)=(x,y)
end


ProductDomain(A,B)=ProductDomain((A,B))
*(A::Domain,B::Domain)=ProductDomain(A,B)

Base.length(d::ProductDomain)=length(d.domains)
Base.transpose(d::ProductDomain)=ProductDomain(d[2],d[1])
Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]
==(d1::ProductDomain,d2::ProductDomain)=d1.domains==d2.domains

Base.first(d::ProductDomain)=(first(d[1]),first(d[2]))

function checkpoints(d::ProductDomain)
    ptsx=checkpoints(d[1])
    ptsy=checkpoints(d[2])
    ret=Array(@compat(Tuple{eltype(d[1]),eltype(d[2])}),0)
    for x in ptsx,y in ptsy
        push!(ret,(x,y))
    end
    ret
end
