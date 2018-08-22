## A domain representing a union of subdomains

export UnionDomain


"""
    UnionDomain((d1,d2,…,dn))

represents a union of multiple subdomains: `{x : x ∈ d1 || … || x ∈ dn}`.
"""
struct UnionDomain{DD,T} <: Domain{T}
    domains::DD
end

UnionDomain(d::Tuple{}) = error("Cannot create UnionDomain with no components")
UnionDomain(d::Tuple) =
    UnionDomain{typeof(d),mapreduce(eltype,promote_type,d)}(d)
UnionDomain(d::AbstractVector) = UnionDomain{typeof(d),eltype(eltype(d))}(d)


UnionDomain(d1::UnionDomain,d2::UnionDomain) = UnionDomain((d1.domains...,d2.domains...))
UnionDomain(d1::UnionDomain{<:AbstractVector},d2::UnionDomain{<:AbstractVector}) =
    UnionDomain([d1.domains ; d2.domains])

UnionDomain(d1::Domain,d2::UnionDomain) = UnionDomain((d1,d2.domains...))
UnionDomain(d1::UnionDomain,d2::Domain) = UnionDomain((d1.domains...,d2))
UnionDomain(d1::Domain,d2::UnionDomain{<:AbstractVector}) = UnionDomain([d1;d2.domains])
UnionDomain(d1::UnionDomain{<:AbstractVector},d2::Domain) = UnionDomain([d1.domains;d2])



UnionDomain(d1::Domain,d2::Domain) = UnionDomain((d1,d2))



canonicaldomain(d::UnionDomain) = d  # we could map all to canonical, but then there would be overlap

isambiguous(d::UnionDomain) = isempty(d.domains)
convert(::Type{UnionDomain{DD,T}},::AnyDomain) where {DD,T} =
    UnionDomain{DD,T}(map(D->D(AnyDomain()),DD.parameters))
convert(::Type{IT},::AnyDomain) where {IT<:UnionDomain} = UnionDomain(tuple())



union(d::Domain) = d
function union(d::AbstractVector{D}) where D<:Domain
    isempty(d) && return EmptyDomain()
    length(d)==1 && return d[1]
    UnionDomain(d)
end
#TODO Check for intersection


union(d1::EmptyDomain,d2::EmptyDomain) = d1
union(d1::EmptyDomain,d2::Domain) = d2
union(d1::Domain,d2::EmptyDomain) = d1

union(d1::AnyDomain,d2::AnyDomain) = d1
union(d1::AnyDomain,d2::Domain) = d2
union(d1::Domain,d2::AnyDomain) = d1

function union(d1::Domain,d2::Domain)
    if d1==d2
        return d1
    end

    Γ=d1∩d2
    if isempty(Γ)
        UnionDomain(d1,d2)
    else
        (d1\Γ) ∪ Γ ∪ (d2\Γ)
    end
end


intersect(d1::UnionDomain,d2::UnionDomain) = mapreduce(d->d1∩d,∪,d2.domains)
intersect(d1::Domain,d2::UnionDomain) = mapreduce(d->d1∩d,∪,d2.domains)
intersect(d1::UnionDomain,d2::Domain) = mapreduce(d->d2∩d,∪,d1.domains)


setdiff(a::UnionDomain,b::UnionDomain) = mapreduce(d->setdiff(d,b),∪,a.domains)
setdiff(a::UnionDomain,b::Domain) = mapreduce(d->setdiff(d,b),∪,a.domains)
setdiff(a::Domain,b::UnionDomain) = mapreduce(d->setdiff(a,d),∩,b.domains)
setdiff(a::UnionDomain,b) = mapreduce(d->setdiff(d,b),∪,a.domains)
setdiff(a,b::UnionDomain) = mapreduce(d->setdiff(a,d),∩,b.domains)

sort(d::UnionDomain;opts...) = UnionDomain(sort(collect(d.domains);opts...))


for op in (:(first),:(last))
    @eval $op(d::UnionDomain) = $op($op(d.domains))
end

#support tuple set
components(d::UnionDomain) = collect(d.domains)
component(d::UnionDomain,k) = d.domains[k]
ncomponents(d::UnionDomain) = length(d.domains)

pieces(d::UnionDomain) = collect(d.domains)
piece(d::UnionDomain,k) = d.domains[k]
npieces(d::UnionDomain) = length(d.domains)


arclength(d::UnionDomain) = mapreduce(arclength,+,d.domains)

==(d1::UnionDomain,d2::UnionDomain) =
    ncomponents(d1) == ncomponents(d2) &&
        all(Bool[component(d1,k) == component(d2,k) for k=1:ncomponents(d1)])


in(x,d::UnionDomain) = any(a->x∈a,d.domains)
issubset(a::Domain,d::UnionDomain) = (a∪d) == d
reverse(d::UnionDomain) = UnionDomain(reverse(map(reverse,d.domains)))

∂(d::UnionDomain) = mapreduce(∂,union,d.domains)

# determine the number of points per piece
function components_npoints(d, n::Int)
    N = ncomponents(d)
    k = n ÷ N
    r = n - N*k
    [fill(k+1, r); fill(k, N-r)]
end


pieces_npoints(d, n::Int) = components_npoints(d, n)

points(d::UnionDomain,n) = vcat(points.(pieces(d), pieces_npoints(d,n))...)

rand(d::UnionDomain) = rand(component(d,rand(1:length(d))))
checkpoints(d::UnionDomain) = mapreduce(checkpoints,union,d.domains)

function merge(d1::UnionDomain, m::Segment)
    ret=d1.domains

    for k=length(ret):-1:1
        it=intersect(ret[k],m)
        if !isempty(it)
            sa=setdiff(ret[k],it)
            m=setdiff(m,it)
            if isempty(sa)
                ret = [ret[1:k-1]...; it; ret[k+1:end]...]
            else
                ret = [ret[1:k-1]...; sa; it; ret[k+1:end]...]
            end
            if isempty(m)
                break
            end
        end
    end
    if !isempty(m)
        ret = [ret...; m]
    end

    UnionDomain(sort!(ret,by=first))
end

function merge(d1::UnionDomain,d2::UnionDomain)
    ret=d1
    for m in d2.domains
        ret=merge(ret,m)
    end
    ret
end


for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::UnionDomain) = UnionDomain(map(a->$op(c,a),d.domains))
        $op(d::UnionDomain,c::Number) = UnionDomain(map(a->$op(a,c),d.domains))
    end
end
/(d::UnionDomain,c::Number) = UnionDomain(map(a->a/c,d.domains))
