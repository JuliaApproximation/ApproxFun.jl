## A domain representing a union of subdomains

canonicaldomain(d::UnionDomain) = d  # we could map all to canonical, but then there would be overlap

isambiguous(d::UnionDomain) = isempty(d.domains)
convert(::Type{UnionDomain{DD,T}},::AnyDomain) where {DD,T} =
    UnionDomain{DD,T}(map(D->D(AnyDomain()),DD.parameters))
convert(::Type{IT},::AnyDomain) where {IT<:UnionDomain} = UnionDomain(tuple())


#support tuple set
components(d::AbstractVector) = d
components(d::UnionDomain) = elements(d)
component(d::UnionDomain,k) = components(d)[k]
ncomponents(d::UnionDomain) = length(elements(d))

pieces(d::UnionDomain) = elements(d)
piece(d::UnionDomain,k) = pieces(d)[k]
npieces(d::UnionDomain) = length(elements(d))


union(::AnyDomain, d::UnionDomain) = d
union(d::UnionDomain, ::AnyDomain) = d

arclength(d::UnionDomain) = mapreduce(arclength,+,d.domains)


reverseorientation(d::UnionDomain) = UnionDomain(reverse(map(reverseorientation,d.domains)))

leftendpoint(d::UnionDomain) = leftendpoint(first(elements(d)))
rightendpoint(d::UnionDomain) = rightendpoint(last(elements(d)))

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


## to move over
function Base.merge(d1::UnionDomain,m::IntervalOrSegment)
    T = float(promote_type(eltype(d1), eltype(m)))
    ret=d1.domains

    for k=length(ret):-1:1
        it=intersect(ret[k],m)
        if arclength(it) ≠ 0
            sa=setdiff(ret[k],it)
            m=setdiff(m,it)
            if arclength(sa) ≤ eps(T)
                ret = [ret[1:k-1]...; it; ret[k+1:end]...]
            else
                ret = [ret[1:k-1]...; sa; it; ret[k+1:end]...]
            end
            if arclength(m) == 0
                break
            end
        end
    end
    if !isempty(m)
        ret = [ret...; m]
    end

    UnionDomain(sort!(ret,by=leftendpoint))
end

function merge(d1::UnionDomain,d2::UnionDomain)
    ret=d1
    for m in d2.domains
        ret=merge(ret,m)
    end
    ret
end
