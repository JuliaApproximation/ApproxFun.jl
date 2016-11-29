export DiracDelta

for TYP in (:DiracSpace,:PointSpace)
    @eval begin
        immutable $TYP{T}<:RealUnivariateSpace{AnyDomain}
          points::Vector{T}
          $TYP(pts::Vector{T})=new(sort(pts))
        end

        $TYP(points::AbstractVector) = $TYP{eltype(points)}(points)
        $TYP() = $TYP(Float64[])
        $TYP(point::Number) = $TYP([point])
        $TYP(p::Point)=$TYP(p.x)
        dimension(d::$TYP)=length(d.points)

        domain(DS::$TYP)=mapreduce(Point,union,DS.points)
        setdomain(DS::$TYP,d::UnionDomain)=$TYP(map(d->d.x,d.domains))
        points(sp::$TYP,n::Integer)=sp.points[1:n]

        spacescompatible(a::$TYP,b::$TYP)=a.points==b.points
        canonicalspace(a::$TYP)=a

        union_rule(a::$TYP,b::$TYP)=$TYP(sort(union(a.points,b.points)))

        function coefficients(cfs::Vector,fromspace::$TYP,tospace::$TYP)
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
                error("The domain of the space being converted from has points that cannot be represented in the space you are converting to.")
            end
        end
    end
end

Space(d::Point)=PointSpace(d)

identity_fun(S::PointSpace)=Fun(S,S.points)
identity_fun(S::DiracSpace)=Fun(PointSpace(S.points),S.points)
transform(S::PointSpace,v::Vector,plan...)=v


function evaluate(f::AbstractVector,PS::PointSpace,x::Number)
    p = findfirst(y->isapprox(x,y),PS.points)
    if p==0
        zero(eltype(f))
    else
        f[p]
    end
end

evaluate(f::AbstractVector,PS::PointSpace,x::AbstractVector)=
    map(y->evaluate(f,PS,y),x)

Base.sum{DS<:DiracSpace}(f::Fun{DS})=sum(f.coefficients[1:dimension(space(f))])



DiracDelta(x::Number)=Fun(DiracSpace(x),[1.])
DiracDelta()=DiracDelta(0.)


function Base.cumsum{S<:DiracSpace,T<:Real}(f::Fun{S},d::Segment{T})
    pts=space(f).points
    @assert pts ==sort(pts)
    cfs=cumsum(f.coefficients)
    if first(d) < first(pts) && last(d) > last(pts)
        Fun(HeavisideSpace([first(d);pts;last(d)]),[0;cfs])
    elseif first(d) == first(pts) && last(d) > last(pts)
        Fun(HeavisideSpace([pts;last(d)]),cfs)
    elseif first(d) < first(pts) && last(d) == last(pts)
        Fun(HeavisideSpace([first(d);pts]),[0;cfs])
    elseif first(d) == first(pts) && last(d) == last(pts)
        Fun(HeavisideSpace(pts),cfs[1:end-1])
    else
        error("Implement")
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


# function Multiplication{PS<:PointSpace}(f::Fun{PS},PS2::PointSpace)
#     @assert space(f).points==PS2.points
#     FiniteOperator(diagm(f.coefficients),PS2,PS2)
# end

# function Multiplication{PS<:PointSpace}(f::Fun{PS},DS::DiracSpace)
#     @assert space(f).points==DS.points
#     FiniteOperator(diagm(f.coefficients),DS,DS)
# end

# function Multiplication{DS<:DiracSpace}(f::Fun{DS},PS::PointSpace)
#     @assert space(f).points==PS.points
#     FiniteOperator(diagm(f.coefficients),PS,space(f))
# end

function coefficienttimes{PS<:PointSpace,DS<:DiracSpace}(f::Fun{PS},g::Fun{DS})
    @assert space(f).points==space(g).points
    Fun(space(g),f.coefficients.*g.coefficients)
end

function coefficienttimes{PS<:PointSpace,DS<:DiracSpace}(f::Fun{DS},g::Fun{PS})
    @assert space(f).points==space(g).points
    Fun(space(f),f.coefficients.*g.coefficients)
end

function coefficienttimes{PS<:PointSpace,PS2<:PointSpace}(f::Fun{PS},g::Fun{PS2})
    @assert space(f).points==space(g).points
    Fun(space(g),f.coefficients.*g.coefficients)
end

./{PS<:PointSpace}(f::Fun,g::Fun{PS}) = f.*inv(g)
Base.inv{PS<:PointSpace}(f::Fun{PS}) = Fun(space(f),1./f.coefficients)
