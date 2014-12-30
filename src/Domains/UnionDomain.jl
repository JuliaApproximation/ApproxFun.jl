## A domain representing a union of subdomains

immutable UnionDomain{D<:Domain} <:Domain
    domains::Vector{D}
end

∪(d1::UnionDomain,d2::UnionDomain)=UnionDomain([d1.domains,d2.domains])
∪(d1::Domain,d2::UnionDomain)=UnionDomain([d1,d2.domains])
∪(d1::UnionDomain,d2::Domain)=UnionDomain([d1.domains,d2])
∪(d1::Domain,d2::Domain)=UnionDomain([d1,d2])
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

    [vcat([points(d.domains[j],k+1) for j=1:r]...),
        vcat([points(d.domains[j],k) for j=r+1:length(d)]...)]
end

Base.rand(d::UnionDomain)=rand(d[rand(1:length(d))])

function Base.merge(d1::UnionDomain,m::Interval)
    ret=d1.domains

    for k=length(ret):-1:1
        it=intersect(ret[k],m)
        if it != []
            sa=setdiff(ret[k],it)
            m=setdiff(m,it)        
            ret = [ret[1:k-1]...,sa,it,ret[k+1:end]...]        
            if m==[]
                break
            end
        end
    end
    @assert m==[]
    UnionDomain(sort!([m,ret...],by=first))
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

