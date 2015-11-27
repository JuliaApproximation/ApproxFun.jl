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

Base.in(x,d::Point)=isapprox(x,d.x)

Base.issubset(a::Point,d::UnionDomain)=a.x in d
Base.issubset(a::Point,b::Domain)=a.x in b

Base.intersect(a::Point,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(a::UnionDomain,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(a::Domain,b::Point)=b.x in a?b:EmptyDomain()
Base.intersect(b::Point,a::UnionDomain)=b.x in a?b:EmptyDomain()
Base.intersect(b::Point,a::Domain)=b.x in a?b:EmptyDomain()

Base.setdiff(a::Point,b::Point)=a==b?EmptyDomain():a
