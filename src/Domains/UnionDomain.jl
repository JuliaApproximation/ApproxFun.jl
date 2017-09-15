## A domain representing a union of subdomains


canonicaldomain(d::UnionDomain) = d  # we could map all to canonical, but then there would be overlap

isambiguous(d::UnionDomain) = isempty(d.domains)
convert(::Type{UnionDomain{DD,T}},::AnyDomain) where {DD,T} =
    UnionDomain{DD,T}(map(D->D(AnyDomain()),DD.parameters))
convert(::Type{IT},::AnyDomain) where {IT<:UnionDomain} = UnionDomain(tuple())

#support tuple set
components(d::UnionDomain) = collect(d.domains)
component(d::UnionDomain,k) = d.domains[k]
ncomponents(d::UnionDomain) = length(d.domains)

pieces(d::UnionDomain) = collect(d.domains)
piece(d::UnionDomain,k) = d.domains[k]
npieces(d::UnionDomain) = length(d.domains)


arclength(d::UnionDomain) = mapreduce(arclength,+,d.domains)

Base.reverse(d::UnionDomain) = UnionDomain(reverse(map(reverse,d.domains)))

# determine the number of points per piece
function pieces_npoints(d, n::Int)
    N = npieces(d)
    k = n ÷ N
    r = n - N*k
    [fill(k+1, r); fill(k, N-r)]
end


points(d::UnionDomain,n) = vcat(points.(pieces(d), pieces_npoints(d,n))...)

Base.rand(d::UnionDomain) = rand(component(d,rand(1:length(d))))
checkpoints(d::UnionDomain) = mapreduce(checkpoints,union,d.domains)
