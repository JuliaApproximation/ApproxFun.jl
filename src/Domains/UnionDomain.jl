## A domain representing a union of subdomains

export UnionDomain


"""
    UnionDomain represents a union of multiple subdomains.
"""
immutable UnionDomain{DD,T,d} <: Domain{T,d}
    domains::DD
end

UnionDomain(d::Tuple)=UnionDomain{typeof(d),mapreduce(eltype,promote_type,d),mapreduce(ndims,max,d)}(d)
UnionDomain(d::Vector)=UnionDomain(tuple(d...))


isambiguous(d::UnionDomain)=isempty(d.domains)
Base.convert{DD,T,d}(::Type{UnionDomain{DD,T,d}},::AnyDomain)=UnionDomain{DD,T,d}(map(D->D(AnyDomain()),DD.parameters))
Base.convert{IT<:UnionDomain}(::Type{IT},::AnyDomain)=UnionDomain(tuple())



∪(d::Domain) = d
∪{D<:Domain}(d::Vector{D}) = UnionDomain(d)
function ∪{D<:Domain}(::Type{D},x)
    out = map(D,x)
    length(out) > 1 ? ∪(out) : out[1]
end
function ∪{D<:Domain}(::Type{D},x,y)
    out = map(D,x,y)
    length(out) > 1 ? ∪(out) : out[1]
end
∪(d1::UnionDomain,d2::UnionDomain)=UnionDomain((d1.domains...,d2.domains...))
∪(d1::Domain,d2::UnionDomain)=UnionDomain((d1,d2.domains...))
∪(d1::UnionDomain,d2::Domain)=UnionDomain((d1.domains...,d2))
∪(d1::Domain,d2::Domain)=UnionDomain((d1,d2))
Base.length(d::UnionDomain)=d.domains|>length
Base.getindex(d::UnionDomain,k)=d.domains[k]
for op in (:(Base.first),:(Base.last))
    @eval $op(d::UnionDomain)=$op($op(d.domains))
end

==(d1::UnionDomain,d2::UnionDomain)=length(d1)==length(d2)&&all(Bool[d1[k]==d2[k] for k=1:length(d1)])


Base.in(x,d::UnionDomain)=any(a->x∈a,d.domains)
Base.issubset(a::Domain,d::UnionDomain)=(a∪d)==d

∂(d::UnionDomain)=mapreduce(∂,union,d.domains)

function points(d::UnionDomain,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d.domains[j],k+1) for j=1:r]...);
        vcat([points(d.domains[j],k) for j=r+1:length(d)]...)]
end

Base.rand(d::UnionDomain)=rand(d[rand(1:length(d))])
checkpoints(d::UnionDomain)=mapreduce(checkpoints,union,d.domains)

function Base.merge(d1::UnionDomain,m::Interval)
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


for op in (:*,:+,:-,:.*,:.+,:.-)
    @eval begin
        $op(c::Number,d::UnionDomain)=UnionDomain(map(a->$op(c,a),d.domains))
        $op(d::UnionDomain,c::Number)=UnionDomain(map(a->$op(a,c),d.domains))
    end
end
/(d::UnionDomain,c::Number) = UnionDomain(map(a->a/c,d.domains))
