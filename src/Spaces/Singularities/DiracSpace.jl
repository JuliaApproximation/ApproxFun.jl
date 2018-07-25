export DiracDelta, KroneckerDelta

for TYP in (:DiracSpace,:PointSpace)
    @eval begin
        struct $TYP{T,D,R} <: Space{D,R}
          points::Vector{T}
          $TYP{T,D,R}(pts::AbstractVector{T}) where {T,D,R} = new(sort(pts))
        end

        function $TYP(points::AbstractVector{T}) where T
            $TYP{eltype(points),UnionDomain{Vector{Point{T}},Point{T}},real(T)}(points)
        end

        $TYP(points::Tuple) = $TYP(collect(points))
        $TYP() = $TYP(Float64[])
        $TYP(point::Number) = $TYP([point])
        $TYP(p::Point)=$TYP(p.x)


        dimension(d::$TYP)=length(d.points)

        # all points are equal, so only one block
        blocklengths(C::$TYP) = [length(C.points)]
        block(C::$TYP,k) = Block(1)

        domain(DS::$TYP) = UnionDomain(Point.(DS.points))
        setdomain(DS::$TYP,d::UnionDomain) = $TYP(map(d->d.x,components(d)))
        points(sp::$TYP,n::Integer)=sp.points[1:n]

        spacescompatible(a::$TYP,b::$TYP) = a.points == b.points
        canonicalspace(a::$TYP) = a

        union_rule(a::$TYP,b::$TYP) = $TYP(sort(union(a.points,b.points)))

        function coefficients(cfs::AbstractVector,fromspace::$TYP,tospace::$TYP)
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

Space(d::Point) = PointSpace(d)

Fun(::typeof(identity), S::PointSpace) = Fun(S,S.points)
Fun(::typeof(identity), S::DiracSpace) = Fun(PointSpace(S.points),S.points)
transform(S::PointSpace,v::AbstractVector,plan...) = v
values(f::Fun{S}) where S<:PointSpace = coefficient(f,:)

function evaluate(f::AbstractVector,PS::PointSpace,x::Number)
    p = findfirst(y->isapprox(x,y),PS.points)
    if p == 0 || p > length(f)
        zero(eltype(f))
    else
        f[p]
    end
end

Base.sum(f::Fun{DS}) where DS<:DiracSpace =
    sum(f.coefficients[1:dimension(space(f))])

DiracDelta(x::Number)=Fun(DiracSpace(x),[1.])
DiracDelta()=DiracDelta(0.)

KroneckerDelta(x::Number) = Fun(PointSpace(x),[1.])
KroneckerDelta() = KroneckerDelta(0.)

function Base.cumsum(f::Fun{S},d::Segment{T}) where {S<:DiracSpace,T<:Real}
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
#     function coefficients(cfs::AbstractVector,fromspace::$TYP,tospace::DiracSpace)
#       [0*tospace.points;coefficients(cfs,fromspace,tospace.space)]
#     end
#   end
# end

# function coefficients(cfs::AbstractVector,fromspace::DiracSpace,tospace::Space)
#   n = length(fromspace.points)
#   if n == 0 || cfs[1:n] == 0*cfs[1:n]
#       coefficients(fromspace.space,tospace)
#   else
#     error("The space you are converting from has Dirac deltas that cannot be represented in the space you are converting to.")
#   end
# end


function Multiplication(f::Fun{PS},PS2::PointSpace) where PS<:PointSpace
    @assert space(f).points==PS2.points
    FiniteOperator(diagm(values(f)),PS2,PS2)
end

function Multiplication(f::Fun{PS},DS::DiracSpace) where PS<:PointSpace
    @assert space(f).points==DS.points
    FiniteOperator(diagm(values(f)),DS,DS)
end

function Multiplication(f::Fun{DS},PS::PointSpace) where DS<:DiracSpace
    @assert space(f).points==PS.points
    FiniteOperator(diagm(coefficient(f,:)),PS,space(f))
end

function coefficienttimes(f::Fun{PS},g::Fun{DS}) where {PS<:PointSpace,DS<:DiracSpace}
    @assert space(f).points==space(g).points
    Fun(space(g),f.coefficients.*g.coefficients)
end

function coefficienttimes(f::Fun{DS},g::Fun{PS}) where {PS<:PointSpace,DS<:DiracSpace}
    @assert space(f).points==space(g).points
    Fun(space(f),f.coefficients.*g.coefficients)
end

function coefficienttimes(f::Fun{PS},g::Fun{PS2}) where {PS<:PointSpace,PS2<:PointSpace}
    @assert space(f).points==space(g).points
    Fun(space(g),f.coefficients.*g.coefficients)
end

/(f::Fun,g::Fun{PS}) where PS<:PointSpace = f*inv(g)
Base.inv(f::Fun{PS}) where PS<:PointSpace = Fun(space(f),1./f.coefficients)

#DiracSpace sampling
function randweights(cfs)
    cfs = f.coefficients/sum(f.coefficients)
    cs = cumsum(cfs)
    r = rand()
    if r≤cs[1]
        return f.space.points[1]
    else
        for n=1:length(cfs)-1
            if cs[n]<r≤cs[n+1] && return f.space.points[n+1]
                break
            end
        end
    end
end


function sample(f::Fun{<:DiracSpace})
    cfs = f.coefficients/sum(f.coefficients)
    randweights(cfs)
end


#integrate DiracSpace
function integrate(f::Fun{<:DiracSpace})
    pts=f.space.points
    cfs=f.coefficients
    int=Fun(HeavisideSpace([pts;Inf]),cumsum(cfs))
    return int
end

