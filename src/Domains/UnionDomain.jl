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




## to move over

Base.setdiff(a::UnionDomain,b::UnionDomain) = mapreduce(d->setdiff(d,b),∪,a.domains)
Base.setdiff(a::UnionDomain,b::Domain) = mapreduce(d->setdiff(d,b),∪,a.domains)
Base.setdiff(a::Domain,b::UnionDomain) = mapreduce(d->setdiff(a,d),∩,b.domains)
Base.setdiff(a::UnionDomain,b) = mapreduce(d->setdiff(d,b),∪,a.domains)
Base.setdiff(a,b::UnionDomain) = mapreduce(d->setdiff(a,d),∩,b.domains)


function Base.merge(d1::UnionDomain,m::Segment)
    ret=d1.domains

    for k=length(ret):-1:1
        it=intersect(ret[k],m)
        if !isempty(it)
            sa=setdiff(ret[k],it)
            m=setdiff(m,it)
            if isempty(sa)
                ret = [ret[1:k-1]...;it;ret[k+1:end]...]
            else
                ret = [ret[1:k-1]...;sa;it;ret[k+1:end]...]
            end
            if isempty(m)
                break
            end
        end
    end
    @assert isempty(m)
    UnionDomain(sort!(ret,by=first))
end

function Base.merge(d1::UnionDomain,d2::UnionDomain)
    ret=d1
    for m in d2.domains
        ret=merge(ret,m)
    end
    ret
end
