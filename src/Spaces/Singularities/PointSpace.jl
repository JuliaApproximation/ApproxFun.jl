immutable PointSpace{T}<:RealUnivariateSpace{AnyDomain}
  points::Vector{T}
  PointSpace(pts::Vector{T})=new(sort(pts))
end

PointSpace(points::AbstractVector) = PointSpace{eltype(points)}(points)
PointSpace() = PointSpace(Float64[])
PointSpace(point::Number) = PointSpace([point])

dimension(d::PointSpace)=length(d.points)

domain(PS::PointSpace)=mapreduce(Point,union,DS.points)
setdomain(PS::PointSpace,d::UnionDomain)=PointSpace(map(d->d.x,d.domains))

spacescompatible(a::PointSpace,b::PointSpace)=a.points==b.points
spacescompatible(a::DiracSpace,b::PointSpace)=a.points==b.points
spacescompatible(a::PointSpace,b::DiracSpace)=a.points==b.points
canonicalspace(a::PointSpace)=a


Base.sum{PS<:PointSpace}(f::Fun{PS})=sum(f.coefficients[1:dimension(space(f))])

union_rule(a::PointSpace,b::PointSpace)=PointSpace(sort(union(a.points,b.points)))

function coefficients(cfs::Vector,fromspace::PointSpace,tospace::PointSpace)
    if spacescompatible(fromspace,tospace)
        return cfs
    end

    @assert length(cfs) ≤ length(fromspace.points)

    # this first for-loop removes coefficients of Dirac points that are zero
    nonzerofromspacepoints = eltype(fromspace.points)[]
    nonzerocfs = eltype(cfs)[]
    for i = 1:length(cfs)
        if cfs[i] != 0
            push!(nonzerofromspacepoints, fromspace.points[i])
            push!(nonzerocfs, cfs[i])
        end
    end


    # if the points that remain can be represented in the tospace
    if issubset(nonzerofromspacepoints,tospace.points)
        finalcfs = eltype(cfs)[]
        j=1 #counter for the nonzerofromspacepoints
        for i = 1:length(tospace.points)
            if j > length(nonzerofromspacepoints)
                break
            elseif nonzerofromspacepoints[j]==tospace.points[i]
                push!(finalcfs,nonzerocfs[j])
                j += 1
            else
                push!(finalcfs,0)
            end
        end
        finalcfs
    else
        error("The point space you are converting from has points that cannot be represented in the point space you are converting to.")
    end
end

function Base.cumsum{S<:PointSpace,T<:Real}(f::Fun{S},d::Interval{T})
  Fun([0.],f.space)
end

function evaluate(f::AbstractVector,PS::PointSpace,x...)
  p = findfirst(x.==PS.points)
  if p == []
    0
  else
    (f.coefficients)[p]
  end
end

# function evaluate(f::AbstractVector,PS::PointSpace,x::Array)
#   output = zeros(size(x))
#   for p in x
#     p = findfirst(x.==PS.points)
#     if p == []
#       output[getindex(p,x)]
#     else
#       (f.coefficients)[p]
#     end
#   end
# end

# function *(f::Fun{DiracSpace},g::Fun{PointSpace})

# end
