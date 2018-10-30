
canonicaldomain(d::ProductDomain) = ProductDomain(map(canonicaldomain,d.domains)...)


# product domains are their own canonical domain
for OP in (:fromcanonical,:tocanonical)
    @eval begin
        $OP(d::ProductDomain, x::Vec) = Vec(map($OP,d.domains,x)...)
        $OP(d::ProductDomain, x::Vec{2}) = Vec($OP(first(d.domains), first(x)), $OP(last(d.domains), last(x)))
    end
end

nfactors(d::ProductDomain) = length(d.domains)
factors(d::ProductDomain) = d.domains
factor(d::ProductDomain,k::Integer) = d.domains[k]

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
    ret=Vector{Vec{length(d.domains),float(mapreduce(eltype,promote_type,d.domains))}}(undef, 0)

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

reverseorientation(d::ProductDomain) = ProductDomain(map(reverseorientation, d.domains))

domainscompatible(a::ProductDomain,b::ProductDomain) =
                        length(a.domains)==length(b.domains) &&
                        all(map(domainscompatible,a.domains,b.domains))
