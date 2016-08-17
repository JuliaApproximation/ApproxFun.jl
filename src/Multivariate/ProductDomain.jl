
export ProductDomain


immutable ProductDomain{D,T,dim} <: Domain{T,dim}
    domains::D
end

ProductDomain(d::Tuple) =
    ProductDomain{typeof(d),mapreduce(eltype,promote_type,d),mapreduce(dimension,+,d)}(d)

fromcanonical(d::BivariateDomain,x::Tuple)=fromcanonical(d,x...)
tocanonical(d::BivariateDomain,x::Tuple)=tocanonical(d,x...)

# product domains are their own canonical domain
for OP in (:fromcanonical,:tocanonical)
    @eval $OP(::ProductDomain,x,y)=(x,y)
end


ProductDomain(A,B)=ProductDomain((A,B))
*(A::ProductDomain,B::ProductDomain)=ProductDomain(tuple(A.domains...,B.domains...))
*(A::ProductDomain,B::Domain)=ProductDomain(tuple(A.domains...,B))
*(A::Domain,B::ProductDomain)=ProductDomain(tuple(A,B.domains...))
*(A::Domain,B::Domain)=ProductDomain(A,B)

Base.length(d::ProductDomain)=length(d.domains)
Base.transpose(d::ProductDomain)=ProductDomain(d[2],d[1])
Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]
==(d1::ProductDomain,d2::ProductDomain)=d1.domains==d2.domains

Base.first(d::ProductDomain)=(first(d[1]),first(d[2]))


function pushappendpts!(ret,xx,pts)
    if isempty(pts)
        push!(ret,xx)
    else
        for x in pts[1]
            pushappendpts!(ret,(xx...,x),pts[2:end])
        end
    end
    ret
end

function checkpoints(d::ProductDomain)
    pts=map(checkpoints,d.domains)
    ret=Array(Tuple{map(eltype,d.domains)...},0)
    pushappendpts!(ret,(),pts)
    ret
end

function points(d::ProductDomain,n::Tuple)
    @assert length(d.domains) == length(n)
    pts=map(points,d.domains,n)
    ret=Array(Tuple{map(eltype,d.domains)...},0)
    pushappendpts!(ret,(),pts)
    ret
end

Base.reverse(d::ProductDomain)=ProductDomain(map(reverse,d.domains))

domainscompatible(a::ProductDomain,b::ProductDomain) =
                        length(a.domains)==length(b.domains) &&
                        all(map(domainscompatible,a.domains,b.domains))
