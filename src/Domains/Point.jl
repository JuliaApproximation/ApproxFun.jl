struct Point{T} <: Domain{T}
    x::T

    Point{T}(x::T) where {T} = new(x)
    Point{T}(x::AnyDomain) where {T} = new(T(NaN))
end

Point(x) = Point{typeof(x)}(x)
Point(::AnyDomain) = Point(NaN)

convert(::Type{Point},::AnyDomain) = Point(NaN)
convert(::Type{Point{T}},::AnyDomain) where T = Point{T}(NaN)

convert(::Type{Number}, d::Point) = d.x
convert(::Type{N}, d::Point) where N<:Number = N(d.x)
Number(d::Point) = convert(Number, d.x)

"""
    Point(x)

represents a single point at `x`.
"""
Point(x)


==(a::Point,b::Point) = (isambiguous(a) && isambiguous(b)) || a.x == b.x

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number,d::Point)=Point($op(c,d.x))
        $op(d::Point,c::Number)=Point($op(d.x,c))
    end
end


/(d::Point,c::Number) = Point(d.x/c)

for op in (:+,:-)
    @eval $op(a::Point,b::Point) = Point($op(a.x,b.x))
end


for op in (:*,:+)
    @eval begin
        $op(a::Point,v::Vector)=map(y->$op(a,y),v)
        $op(v::Vector,a::Point)=map(y->$op(y,a),v)
    end
end


isambiguous(d::Point) = isnan(d.x)




norm(p::Point) = norm(p.x)

getindex(p::Point,k...) = p.x[k...]
first(p::Point) = p.x
last(p::Point) = p.x

in(x,d::Point) = isapprox(x,d.x)

isnan(d::Point) = false

issubset(a::Point,d::UnionDomain) = a.x in d
issubset(a::Point,b::Domain) = a.x in b

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


convert(::Type{Domain},a::Number) = Point(a)
