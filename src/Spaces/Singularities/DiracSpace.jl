immutable DiracSpace{T<:Space}<:RealUnivariateSpace{AnyDomain}
  points::Vector{Float64}
end

dimension(d::DiracSpace)=length(d.points)

DiracSpace() = DiracSpace(Float64[])
DiracSpace(points) = DiracSpace(sort(points))

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

    if length(cfs) < length(tospace.points)
        error("You have given fewer coefficients than Dirac points for the space you are converting to.")
    end
    if length(cfs) < length(fromspace.points)
        error("You have given fewer coefficients than Dirac points for the space you are converting from.")
    end
    # this first for-loop removes coefficients of Dirac points that are zero
    nonzerofromspacepoints = []
    nonzeroDiraccfs = []
    for i = 1:length(fromspace.points)
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
            if nonzerofromspacepoints[j] in tospace.points
                push!(finalDiraccfs,nonzeroDirccfs[j])
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
