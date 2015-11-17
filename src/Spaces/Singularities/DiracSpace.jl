export DiracDelta


immutable DiracSpace{T}<:RealUnivariateSpace{AnyDomain}
  points::Vector{T}
  DiracSpace(pts::Vector{T})=new(sort(pts))
end

DiracSpace(points::AbstractVector) = DiracSpace{eltype(points)}(points)
DiracSpace() = DiracSpace(Float64[])
DiracSpace(point::Number) = DiracSpace([point])

dimension(d::DiracSpace)=length(d.points)



#to be extended to include dirac points
domain(DS::DiracSpace)=mapreduce(Point,union,DS.points)
setdomain(DS::DiracSpace,d::UnionDomain)=DiracSpace(map(d->d.x,d.domains))

spacescompatible(a::DiracSpace,b::DiracSpace)=a.points==b.points
canonicalspace(a::DiracSpace)=a

Base.sum{DS<:DiracSpace}(f::Fun{DS})=sum(f.coefficients[1:dimension(space(f))])

union_rule(a::DiracSpace,b::DiracSpace)=DiracSpace(sort(union(a.points,b.points)))

function coefficients(cfs::Vector,fromspace::DiracSpace,tospace::DiracSpace)
    if spacescompatible(fromspace,tospace)
        return cfs
    end

    @assert length(cfs) ≤ length(fromspace.points)

    # this first for-loop removes coefficients of Dirac points that are zero
    nonzerofromspacepoints = eltype(fromspace.points)[]
    nonzeroDiraccfs = eltype(cfs)[]
    for i = 1:length(cfs)
        if cfs[i] != 0
            push!(nonzerofromspacepoints, fromspace.points[i])
            push!(nonzeroDiraccfs, cfs[i])
        end
    end


    # if the points that remain can be represented in the tospace
    if issubset(nonzerofromspacepoints,tospace.points)
        finalDiraccfs = eltype(cfs)[]
        j=1 #counter for the nonzerofromspacepoints
        for i = 1:length(tospace.points)
            if j > length(nonzerofromspacepoints)
                break
            elseif nonzerofromspacepoints[j]==tospace.points[i]
                push!(finalDiraccfs,nonzeroDiraccfs[j])
                j += 1
            else
                push!(finalDiraccfs,0)
            end
        end
        finalDiraccfs
    else
        error("The space you are converting from has Dirac deltas that cannot be represented in the space you are converting to.")
    end
end

DiracDelta(x::Number)=Fun([1.],DiracSpace(x))
DiracDelta()=DiracDelta(0.)

# for TYP in (:ReSpace,:Space)
#   @eval begin
#     function coefficients(cfs::Vector,fromspace::$TYP,tospace::DiracSpace)
#       [0*tospace.points;coefficients(cfs,fromspace,tospace.space)]
#     end
#   end
# end

# function coefficients(cfs::Vector,fromspace::DiracSpace,tospace::Space)
#   n = length(fromspace.points)
#   if n == 0 || cfs[1:n] == 0*cfs[1:n]
#       coefficients(fromspace.space,tospace)
#   else
#     error("The space you are converting from has Dirac deltas that cannot be represented in the space you are converting to.")
#   end
# end
