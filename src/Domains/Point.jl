immutable Point{T} <: Domain{T,0}
    x::T
end

==(a::Point,b::Point)=a.x==b.x

for op in (:*,:+,:-,:.*,:.+,:.-,:.^)
    @eval begin
        $op(c::Number,d::Point)=Point($op(c,d.x))
        $op(d::Point,c::Number)=Point($op(d.x,c))
    end
end

for op in (:/,:./)
    @eval $op(d::Point,c::Number)=Point($op(d.x,c))
end

for op in (:+,:-,:.+,:.-)
    @eval $op(a::Point,b::Point)=Point($op(a.x,b.x))
end


for op in (:*,:+)
    @eval begin
        $op(a::Point,v::Vector)=map(y->$op(a,y),v)
        $op(v::Vector,a::Point)=map(y->$op(y,a),v)
    end
end


isambiguous(d::Point) = isnan(d.x)
Base.convert{PT<:Point}(::Type{PT},::AnyDomain) = PT(NaN)

Base.norm(p::Point)=norm(p.x)

Base.getindex(p::Point,k...)=p.x[k...]

Base.in(x,d::Point) = isapprox(x,d.x)

Base.isnan(d::Point) = false

Base.issubset(a::Point,d::UnionDomain)=a.x in d
Base.issubset(a::Point,b::Domain)=a.x in b

Base.intersect(a::Point,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(a::UnionDomain,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(a::Domain,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(b::Point,a::UnionDomain)=b.x in a?b:EmptyDomain()
Base.intersect(b::Point,a::Domain)=b.x in a?b:EmptyDomain()

Base.setdiff(a::Point,b::Point)=a==b?EmptyDomain():a
Base.reverse(a::Point)=a


canonicaldomain(a::Point)=Point(0.)

points(a::Point,n) = eltype(a)[a.x]
checkpoints(a::Point) = eltype(a)[a.x]
