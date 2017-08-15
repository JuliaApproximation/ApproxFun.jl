


export JacobiWeight



"""
    JacobiWeight(β,α,s::Space)

weights a space `s` by a Jacobi weight, which on `-1..1`
is `(1+x)^β*(1-x)^α`.
For other domains, the weight is inferred by mapping to `-1..1`.
"""
#TODO: support general types for paramaters
struct JacobiWeight{S,DD,RR} <: WeightSpace{S,DD,RR}
    β::Float64
    α::Float64
    space::S
    function JacobiWeight{S,DD,RR}(β::Float64,α::Float64,space::S) where {S,DD,RR}
        if isa(space,JacobiWeight)
            new(β+space.β,α+space.α,space.space) else
            new(β,α,space)
        end
    end
end

JacobiWeight(a::Number,b::Number,d::Space) =
    JacobiWeight{typeof(d),domaintype(d),rangetype(d)}(Float64(a),Float64(b),d)
JacobiWeight(β::Number,α::Number,d::JacobiWeight) =  JacobiWeight(β+d.β,α+d.α,d.space)
JacobiWeight(a::Number,b::Number,d::IntervalDomain) = JacobiWeight(Float64(a),Float64(b),Space(d))
JacobiWeight(a::Number,b::Number,d) = JacobiWeight(Float64(a),Float64(b),Space(d))
JacobiWeight(a::Number,b::Number) = JacobiWeight(a,b,Chebyshev())

JacobiWeight(a::Number,b::Number,s::PiecewiseSpace) = PiecewiseSpace(JacobiWeight(a,b,components(s)))

identity_fun(S::JacobiWeight)=isapproxinteger(S.β)&&isapproxinteger(S.α)?Fun(x->x,S):Fun(identity,domain(S))

order(S::JacobiWeight{Ultraspherical{Int,D,R},D,R}) where {D,R} = order(S.space)


spacescompatible(A::JacobiWeight,B::JacobiWeight) =
    A.β ≈ B.β && A.α ≈ B.α && spacescompatible(A.space,B.space)
spacescompatible(A::JacobiWeight,B::Space{DD,RR}) where {DD<:IntervalDomain,RR<:Real} =
    spacescompatible(A,JacobiWeight(0,0,B))
spacescompatible(B::Space{DD,RR},A::JacobiWeight) where {DD<:IntervalDomain,RR<:Real} =
    spacescompatible(A,JacobiWeight(0,0,B))

transformtimes(f::Fun{JW1},g::Fun{JW2}) where {JW1<:JacobiWeight,JW2<:JacobiWeight}=
            Fun(JacobiWeight(f.space.β+g.space.β,f.space.α+g.space.α,f.space.space),
                coefficients(transformtimes(Fun(f.space.space,f.coefficients),
                                            Fun(g.space.space,g.coefficients))))
transformtimes(f::Fun{JW},g::Fun) where {JW<:JacobiWeight} =
    Fun(f.space,coefficients(transformtimes(Fun(f.space.space,f.coefficients),g)))
transformtimes(f::Fun,g::Fun{JW}) where {JW<:JacobiWeight} =
    Fun(g.space,coefficients(transformtimes(Fun(g.space.space,g.coefficients),f)))

jacobiweight(β,α,x) = (1+x).^β.*(1-x).^α
jacobiweight(β,α,d::Domain) = Fun(JacobiWeight(β,α,ConstantSpace(d)),[1.])
jacobiweight(β,α) = jacobiweight(β,α,Interval())

weight(sp::JacobiWeight,x) = jacobiweight(sp.β,sp.α,tocanonical(sp,x))
dimension(sp::JacobiWeight) = dimension(sp.space)


Base.first(f::Fun{JW}) where {JW<:JacobiWeight} = space(f).β>0?zero(eltype(f)):f(first(domain(f)))
Base.last(f::Fun{JW}) where {JW<:JacobiWeight} = space(f).α>0?zero(eltype(f)):f(last(domain(f)))

setdomain(sp::JacobiWeight,d::Domain)=JacobiWeight(sp.β,sp.α,setdomain(sp.space,d))

# we assume that points avoids singularities


##TODO: paradigm for same space
function coefficients(f::AbstractVector,sp1::JacobiWeight{SJ1,DD},sp2::JacobiWeight{SJ2,DD}) where {SJ1,SJ2,DD<:IntervalDomain}
    β,α=sp1.β,sp1.α
    c,d=sp2.β,sp2.α

    if isapprox(c,β) && isapprox(d,α)
        # remove wrapper spaces and then convert
        coefficients(f,sp1.space,sp2.space)
    else
        # go back to default
        defaultcoefficients(f,sp1,sp2)
    end
end
coefficients(f::AbstractVector,sp::JacobiWeight{SJ,DD},
             S2::SubSpace{S,IT,DD,RR}) where {SJ,S,IT,DD<:IntervalDomain,RR<:Real} = subspace_coefficients(f,sp,S2)
coefficients(f::AbstractVector,S2::SubSpace{S,IT,DD,RR},
             sp::JacobiWeight{SJ,DD}) where {SJ,S,IT,DD<:IntervalDomain,RR<:Real} = subspace_coefficients(f,S2,sp)
#TODO: it could be possible that we want to JacobiWeight a SumSpace....
coefficients(f::AbstractVector,sp::JacobiWeight{SJ,DD},S2::SumSpace{SV,DD,RR}) where {SJ,SV,DD<:IntervalDomain,RR<:Real} =
    sumspacecoefficients(f,sp,S2)
coefficients(f::AbstractVector,sp::JacobiWeight{SJ,Segment{Vec{2,TT}}},S2::TensorSpace{SV,TTT,DD}) where {SJ,TT,SV,TTT,DD<:BivariateDomain} =
    coefficients(f,sp,JacobiWeight(0,0,S2))

coefficients(f::AbstractVector,sp::JacobiWeight{SJ,DD},S2::Space{DD,RR}) where {SJ,DD<:IntervalDomain,RR<:Real} =
    coefficients(f,sp,JacobiWeight(0,0,S2))
coefficients(f::AbstractVector,sp::ConstantSpace{DD},ts::JacobiWeight{SJ,DD}) where {SJ,DD<:IntervalDomain} =
    f.coefficients[1]*ones(ts).coefficients
coefficients(f::AbstractVector,S2::Space{DD,RR},sp::JacobiWeight{SJ,DD}) where {SJ,DD<:IntervalDomain,RR<:Real} =
    coefficients(f,JacobiWeight(0,0,S2),sp)


"""
`increase_jacobi_parameter(f)` multiplies by `1-x^2` on the unit interval.
`increase_jacobi_parameter(-1,f)` multiplies by `1+x` on the unit interval.
`increase_jacobi_parameter(+1,f)` multiplies by `1-x` on the unit interval.
On other domains this is accomplished by mapping to the unit interval.
"""
increase_jacobi_parameter(f) = Fun(f,JacobiWeight(f.space.β+1,f.space.α+1,space(f).space))
increase_jacobi_parameter(s,f) = s==-1?Fun(f,JacobiWeight(f.space.β+1,f.space.α,space(f).space)):
                                       Fun(f,JacobiWeight(f.space.β,f.space.α+1,space(f).space))



function canonicalspace(S::JacobiWeight)
    if isapprox(S.β,0) && isapprox(S.α,0)
        canonicalspace(S.space)
    else
        #TODO: promote singularities?
        JacobiWeight(S.β,S.α,canonicalspace(S.space))
    end
end

function union_rule(A::ConstantSpace,B::JacobiWeight{P}) where P<:PolynomialSpace
    # we can convert to a space that contains contants provided
    # that the parameters are integers
    # when the parameters are -1 we keep them
    if isapproxinteger(B.β) && isapproxinteger(B.α)
        JacobiWeight(min(B.β,0.),min(B.α,0.),B.space)
    else
        NoSpace()
    end
end


## Algebra

function /(c::Number,f::Fun{JW}) where JW<:JacobiWeight
    g=c/Fun(space(f).space,f.coefficients)
    Fun(JacobiWeight(-f.space.β,-f.space.α,space(g)),g.coefficients)
end

function ^(f::Fun{JW},k::Float64) where JW<:JacobiWeight
    S=space(f)
    g=Fun(S.space,coefficients(f))^k
    Fun(JacobiWeight(k*S.β,k*S.α,space(g)),coefficients(g))
end

function *(f::Fun{JW1},g::Fun{JW2}) where {JW1<:JacobiWeight,JW2<:JacobiWeight}
    @assert domainscompatible(f,g)
    fβ,fα=f.space.β,f.space.α
    gβ,gα=g.space.β,g.space.α
    m=(Fun(space(f).space,f.coefficients).*Fun(space(g).space,g.coefficients))
    if isapprox(fβ+gβ,0)&&isapprox(fα+gα,0)
        m
    else
        Fun(JacobiWeight(fβ+gβ,fα+gα,space(m)),m.coefficients)
    end
end


/(f::Fun{JW1},g::Fun{JW2}) where {JW1<:JacobiWeight,JW2<:JacobiWeight}=f*(1/g)

# O(min(m,n)) Ultraspherical conjugated inner product

function conjugatedinnerproduct(sp::Ultraspherical,u::AbstractVector{S},v::AbstractVector{V}) where {S,V}
    λ=order(sp)
    if λ==1
        mn = min(length(u),length(v))
        if mn > 0
            return dotu(u[1:mn],v[1:mn])*π/2
        else
            return zero(promote_type(eltype(u),eltype(v)))
        end
    else
        T,mn = promote_type(S,V),min(length(u),length(v))
        if mn > 1
            wi = sqrt(convert(T,π))*gamma(λ+one(T)/2)/gamma(λ+one(T))
            ret = u[1]*wi*v[1]
            for i=2:mn
              wi *= (i-2one(T)+2λ)/(i-one(T)+λ)*(i-2one(T)+λ)/(i-one(T))
              ret += u[i]*wi*v[i]
            end
            return ret
        elseif mn > 0
            wi = sqrt(convert(T,π))*gamma(λ+one(T)/2)/gamma(λ+one(T))
            return u[1]*wi*v[1]
        else
            return zero(promote_type(eltype(u),eltype(v)))
        end
    end
end

function conjugatedinnerproduct(::Chebyshev,u::AbstractVector,v::AbstractVector)
    mn = min(length(u),length(v))
    if mn > 1
        return (2u[1]*v[1]+dotu(u[2:mn],v[2:mn]))*π/2
    elseif mn > 0
        return u[1]*v[1]*π
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end


function bilinearform(f::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}},g::Fun{Ultraspherical{LT,D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f).space)
    if order(space(g)) == λ && f.space.β == f.space.α == λ-0.5
        return complexlength(d)/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function bilinearform(f::Fun{Ultraspherical{LT,D,R}},
                    g::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f))
    if order(space(g).space) == λ && g.space.β == g.space.α == λ-0.5
        return complexlength(d)/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function bilinearform(f::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}},
                    g::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f).space)
    if order(space(g).space) == λ && f.space.β+g.space.β == f.space.α+g.space.α == λ-0.5
        return complexlength(domain(f))/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function linebilinearform(f::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}},g::Fun{Ultraspherical{LT,D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f).space)
    if order(space(g)) == λ && f.space.β == f.space.α == λ-0.5
        return arclength(d)/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end

function linebilinearform(f::Fun{Ultraspherical{LT,D,R}},g::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f))
    if order(space(g).space) == λ &&  g.space.β == g.space.α == λ-0.5
        return arclength(d)/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end

function linebilinearform(f::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}},g::Fun{JacobiWeight{Ultraspherical{LT,D,R},D,R}}) where {LT,D,R}
    d = domain(f)
    @assert d == domain(g)
    λ = order(space(f).space)
    if order(space(g).space) == λ &&  f.space.β+g.space.β == f.space.α+g.space.α == λ-0.5
        return arclength(d)/2*conjugatedinnerproduct(Ultraspherical(λ,d),f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end
