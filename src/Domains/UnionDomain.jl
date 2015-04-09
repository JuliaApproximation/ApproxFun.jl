## A domain representing a union of subdomains

export UnionDomain

immutable UnionDomain{D<:Domain,T<:Number} <:Domain{T}
    domains::Vector{D}
end

UnionDomain{D<:Domain}(d::Vector{D})=UnionDomain{D,mapreduce(eltype,promote_type,d)}(d)


∪(d::Domain) = d
∪{D<:Domain}(d::Vector{D}) = UnionDomain(d)
∪{D1,D2,T1,T2}(d1::UnionDomain{D1,T1},d2::UnionDomain{D2,T2})=UnionDomain([d1.domains,d2.domains])
∪{T1,D2,T2}(d1::Domain{T1},d2::UnionDomain{D2,T2})=UnionDomain([d1,d2.domains])
∪{D1,T1,T2}(d1::UnionDomain{D1,T1},d2::Domain{T2})=UnionDomain([d1.domains,d2])
∪{T1,T2}(d1::Domain{T1},d2::Domain{T2})=UnionDomain([d1,d2])
Base.length(d::UnionDomain)=d.domains|>length
Base.getindex(d::UnionDomain,k)=d.domains[k]
for op in (:(Base.first),:(Base.last))
    @eval $op(d::UnionDomain)=d.domains|>$op|>$op
end

==(d1::UnionDomain,d2::UnionDomain)=length(d1)==length(d2)&&all(Bool[d1[k]==d2[k] for k=1:length(d1)])


∂(d::UnionDomain)=mapreduce(∂,union,d.domains)

function points(d::UnionDomain,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d.domains[j],k+1) for j=1:r]...);
        vcat([points(d.domains[j],k) for j=r+1:length(d)]...)]
end

Base.rand(d::UnionDomain)=rand(d[rand(1:length(d))])

checkpoints(d::UnionDomain)=mapreduce(checkpoints,union,d.domains)

function Base.merge{D}(d1::UnionDomain{D},m::Interval)
    ret=d1.domains
    T=promote_type(D,typeof(m))

    for k=length(ret):-1:1
        it=intersect(ret[k],m)
        if !isempty(it)
            sa=setdiff(ret[k],it)
            m=setdiff(m,it)
            if isempty(sa)
                ret = T[ret[1:k-1]...;it;ret[k+1:end]...]
            else
                ret = T[ret[1:k-1]...;sa;it;ret[k+1:end]...]
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

