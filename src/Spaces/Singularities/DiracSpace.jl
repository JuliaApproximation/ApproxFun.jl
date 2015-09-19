immutable DiracSpace{T<:Space}<:RealUnivariateSpace
  space::T
  points::Vector{Float64}
end

DiracSpace() = DiracSpace{Chebyshev}(Chebyshev(),[])
DiracSpace(points) = DiracSpace{Chebyshev}(Chebyshev(),sort(points))

#to be extended to include dirac points
domain(DS::DiracSpace)=domain(DS.space)∪mapreduce(Point,union,DS.points)
setdomain(DS::DiracSpace,d::UnionDomain)=DiracSpace(setdomain(DS.space,first(d.domains)),map(d->d.x,d.domains[2:end]))

spacescompatible(a::DiracSpace,b::DiracSpace)=spacescompatible(a.space,b.space) && a.points==b.points
canonicalspace(a::DiracSpace)=a

function Base.sum{DS<:DiracSpace}(f::Fun{DS})
    n = length(space(f).points)
    sum(f.coefficients[1:n])+sum(Fun(f.coefficients[n+1:end],space(f).space))
end

function evaluate{DS<:DiracSpace}(f::Fun{DS},x::Number)
  n = length(f.space.points)
  if x in f.space.points
    error("You cannot evaluate a Dirac delta at its center.")
  else
    evaluate(Fun(f.coefficients[n+1:end],f.space.space),x)
  end
end

evaluate(S::DiracSpace,coeffs::Vector,x::Number)=evaluate(Fun(coeffs,S),x)

function union_rule(a::DiracSpace,b::DiracSpace)
  DiracSpace(union(a.space,b.space),sort(union(a.points,b.points)))
end
function union_rule(a::DiracSpace,b::Space)
  DiracSpace(union(a.space,b),a.points)
end

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
    finalDiraccfs = []
    j=1 #counter for the nonzerofromspacepoints
    for i = 1:length(tospace.points)
      if nonzerofromspacepoints[j] in tospace.points
        finalDiraccfs[end] = nonzeroDirccfs[j]
        j += 1
      else
        finalDiraccfs[end]=0
      end
    end
    [finalDiraccfs;coefficients(cfs[n+1:end],fromspace.space,tospace.space)]
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
