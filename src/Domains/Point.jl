Point(::AnyDomain) = Point(NaN)
Point{T}(x::AnyDomain) where {T} = new(T(NaN))

convert(::Type{Point},::AnyDomain) = Point(NaN)
convert(::Type{Point{T}},::AnyDomain) where T = Point{T}(NaN)



isambiguous(d::Point) = isnan(d.x)




norm(p::Point) = norm(p.x)

getindex(p::Point,k...) = p.x[k...]
first(p::Point) = p.x
last(p::Point) = p.x

Base.isnan(d::Point) = false


issubset(a::Point,d::UnionDomain) = a.x in d

intersect(a::Point,b::Point) = b.x in a ? b : EmptyDomain()
intersect(a::UnionDomain,b::Point) = b.x in a ? b : EmptyDomain()
intersect(a::Domain,b::Point) = b.x in a ? b : EmptyDomain()
intersect(b::Point,a::UnionDomain) = b.x in a ? b : EmptyDomain()
intersect(b::Point,a::Domain) = b.x in a ? b : EmptyDomain()

setdiff(a::Point,b::Point) = a==b ? EmptyDomain() : a
reverse(a::Point) = a


canonicaldomain(a::Point) = Point(0.)
tocanonical(a::Point,x) = x-a.x
fromcanonical(a::Point,x) = x+a.x

points(a::Point,n) = eltype(a)[a.x]
checkpoints(a::Point) = eltype(a)[a.x]
