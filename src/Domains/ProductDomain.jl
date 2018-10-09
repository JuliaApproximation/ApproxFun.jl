
export ProductDomain

"""
    ProductDomain((d1,d2))

represents the product of two domains, the set `{(x,y) : x ∈ d1 & y ∈ d2}`.

Multiplication of domains is overrident to return a `ProductDomain`.
For example, the following represents the rectangle `1 ≤ x ≤ 2 & 3 ≤ y ≤ 4`:
```julia
Interval(1,2)*(3,4)
```
"""
struct ProductDomain{D,T} <: Domain{T}
    domains::D
end

ProductDomain(d::Tuple) =
    ProductDomain{typeof(d),Vec{length(d),mapreduce(eltype,promote_type,d)}}(d)

issubset(a::ProductDomain,b::ProductDomain) =
  length(a.domains) == length(b.domains) && all(issubset(a.domains[i],b.domains[i]) for i in eachindex(a.domains))


canonicaldomain(d::ProductDomain) = ProductDomain(map(canonicaldomain,d.domains))

# product domains are their own canonical domain
for OP in (:fromcanonical,:tocanonical)
    @eval begin
        $OP(d::ProductDomain, x::Vec) = Vec(map($OP,d.domains,x)...)
        $OP(d::ProductDomain, x::Vec{2}) = Vec($OP(first(d.domains), first(x)), $OP(last(d.domains), last(x)))
    end
end


ProductDomain(A,B) = ProductDomain((A,B))
*(A::ProductDomain, B::ProductDomain) = ProductDomain(tuple(A.domains...,B.domains...))
*(A::ProductDomain, B::Domain) = ProductDomain(tuple(A.domains...,B))
*(A::Domain, B::ProductDomain) = ProductDomain(tuple(A,B.domains...))
*(A::Domain, B::Domain) = ProductDomain(A,B)
^(A::Domain, p::Integer) = p == 1 ? A : A*A^(p-1)


transpose(d::ProductDomain) = ProductDomain(d[2],d[1])
nfactors(d::ProductDomain) = length(d.domains)
factor(d::ProductDomain, k::Integer) = d.domains[k]
==(d1::ProductDomain, d2::ProductDomain) = d1.domains==d2.domains

first(d::ProductDomain) = (first(d[1]),cfirst(d[2]))

in(x::Vec, d::ProductDomain) = reduce(&,map(in,x,d.domains))


function pushappendpts!(ret, xx, pts)
    if isempty(pts)
        push!(ret,Vec(xx...))
    else
        for x in pts[1]
            pushappendpts!(ret,(xx...,x),pts[2:end])
        end
    end
    ret
end

function checkpoints(d::ProductDomain)
    pts=map(checkpoints,d.domains)
    ret=Vector{Vec{length(d.domains),mapreduce(eltype,promote_type,d.domains)}}(undef, 0)

    pushappendpts!(ret,(),pts)
    ret
end

function points(d::ProductDomain,n::Tuple)
    @assert length(d.domains) == length(n)
    pts=map(points,d.domains,n)
    ret=Vector{Vec{length(d.domains),mapreduce(eltype,promote_type,d.domains)}}(undef, 0)
    pushappendpts!(ret,Vec(x),pts)
    ret
end

reverse(d::ProductDomain) = ProductDomain(map(reverse,d.domains))

domainscompatible(a::ProductDomain,b::ProductDomain) =
                        length(a.domains)==length(b.domains) &&
                        all(map(domainscompatible,a.domains,b.domains))
