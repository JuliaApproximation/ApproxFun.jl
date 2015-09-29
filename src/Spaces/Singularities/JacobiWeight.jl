


export JacobiWeight



"""
`JacobiWeight`
weights a basis on `[-1,1]` weighted by `(1+x)^α*(1-x)^β`.
Note the inconsistency of the parameters with `Jacobi`.
when the domain is `[a,b]` the weight is inferred by mapping to `[-1,1]`
"""
immutable JacobiWeight{S,DD} <: WeightSpace{RealBasis,DD,1}
    α::Float64
    β::Float64
    space::S
    function JacobiWeight(α::Float64,β::Float64,space::S)
        if isa(space,JacobiWeight)
            JacobiWeight(α+space.α,β+space.β,space.space)
        else
            new(α,β,space)
        end
    end
end

JacobiWeight(a::Number,b::Number,d::RealUnivariateSpace)=JacobiWeight{typeof(d),typeof(domain(d))}(Float64(a),Float64(b),d)
JacobiWeight(a::Number,b::Number,d::IntervalDomain)=JacobiWeight(Float64(a),Float64(b),Space(d))
JacobiWeight(a::Number,b::Number,d::Vector)=JacobiWeight(Float64(a),Float64(b),Space(d))
JacobiWeight(a::Number,b::Number)=JacobiWeight(a,b,Chebyshev())

JacobiWeight(a::Number,b::Number,s::Vector) = map(s->JacobiWeight(a,b,s),s)
JacobiWeight(a::Number,b::Number,s::PiecewiseSpace) = PiecewiseSpace(JacobiWeight(a,b,vec(s)))

identity_fun(S::JacobiWeight)=isapproxinteger(S.α)&&isapproxinteger(S.β)?Fun(x->x,S):Fun(identity,domain(S))



spacescompatible(A::JacobiWeight,B::JacobiWeight)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)
spacescompatible{DD<:Interval}(A::JacobiWeight,B::RealUnivariateSpace{DD})=spacescompatible(A,JacobiWeight(0,0,B))
spacescompatible{DD<:Interval}(B::RealUnivariateSpace{DD},A::JacobiWeight)=spacescompatible(A,JacobiWeight(0,0,B))

transformtimes{JW1<:JacobiWeight,JW2<:JacobiWeight}(f::Fun{JW1},g::Fun{JW2})=
            Fun(coefficients(transformtimes(Fun(f.coefficients,f.space.space),
                                            Fun(g.coefficients,g.space.space))),
                             JacobiWeight(f.space.α+g.space.α,f.space.β+g.space.β,f.space.space))
transformtimes{JW<:JacobiWeight}(f::Fun{JW},g::Fun) = Fun(coefficients(transformtimes(Fun(f.coefficients,f.space.space),g)),f.space)
transformtimes{JW<:JacobiWeight}(f::Fun,g::Fun{JW}) = Fun(coefficients(transformtimes(Fun(g.coefficients,g.space.space),f)),g.space)

##  α and β are opposite the convention for Jacobi polynomials
# Here, α is the left algebraic singularity and β is the right algebraic singularity.

jacobiweight(α,β,x)=(1+x).^α.*(1-x).^β
weight(sp::JacobiWeight,x)=jacobiweight(sp.α,sp.β,tocanonical(sp,x))



setdomain(sp::JacobiWeight,d::Domain)=JacobiWeight(sp.α,sp.β,setdomain(sp.space,d))

# we assume that points avoids singularities


##TODO: paradigm for same space
function coefficients{SJ1,SJ2,DD<:Interval}(f::Vector,sp1::JacobiWeight{SJ1,DD},sp2::JacobiWeight{SJ2,DD})
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β

    if isapprox(c,α) && isapprox(d,β)
        coefficients(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end
coefficients{SJ,S,n,st,DD<:Interval}(f::Vector,sp::JacobiWeight{SJ,DD},S2::SliceSpace{n,st,S,RealBasis,DD,1})=error("Implement")
coefficients{SJ,S,n,st,DD<:Interval}(f::Vector,S2::SliceSpace{n,st,S,RealBasis,DD,1},sp::JacobiWeight{SJ,DD})=error("Implement")
#TODO: it could be possible that we want to JacobiWeight a SumSpace....
coefficients{SJ,SV,DD<:Interval}(f::Vector,sp::JacobiWeight{SJ,DD},S2::SumSpace{SV,RealBasis,DD,1})=sumspacecoefficients(f,sp,S2)
coefficients{SJ,DD<:Interval}(f::Vector,sp::JacobiWeight{SJ,DD},S2::RealUnivariateSpace{DD})=coefficients(f,sp,JacobiWeight(0,0,S2))
coefficients{SJ,DD<:Interval}(f::Vector,S2::RealUnivariateSpace{DD},sp::JacobiWeight{SJ,DD})=coefficients(f,JacobiWeight(0,0,S2),sp)


"""
`increase_jacobi_parameter(f)` multiplies by `1-x^2` on the unit interval.
`increase_jacobi_parameter(-1,f)` multiplies by `1+x` on the unit interval.
`increase_jacobi_parameter(+1,f)` multiplies by `1-x` on the unit interval.
On other domains this is accomplished by mapping to the unit interval.
"""
increase_jacobi_parameter(f)=Fun(f,JacobiWeight(f.space.α+1,f.space.β+1,space(f).space))
increase_jacobi_parameter(s,f)=s==-1?Fun(f,JacobiWeight(f.space.α+1,f.space.β,space(f).space)):
                                     Fun(f,JacobiWeight(f.space.α,f.space.β+1,space(f).space))



function canonicalspace(S::JacobiWeight)
    if isapprox(S.α,0) && isapprox(S.β,0)
        canonicalspace(S.space)
    else
        #TODO: promote singularities?
        JacobiWeight(S.α,S.β,canonicalspace(S.space))
    end
end

function union_rule{P<:PolynomialSpace}(A::ConstantSpace,B::JacobiWeight{P})
    # we can convert to a space that contains contants provided
    # that the parameters are integers
    # when the parameters are -1 we keep them
    if isapproxinteger(B.α) && isapproxinteger(B.β)
        JacobiWeight(min(B.α,0.),min(B.β,0.),B.space)
    else
        NoSpace()
    end
end


## Algebra

for op in (:/,:./)
    @eval begin
        function ($op){JW<:JacobiWeight}(c::Number,f::Fun{JW})
            g=($op)(c,Fun(f.coefficients,space(f).space))
            Fun(g.coefficients,JacobiWeight(-f.space.α,-f.space.β,space(g)))
        end
    end
end

function .^{JW<:JacobiWeight}(f::Fun{JW},k::Float64)
    S=space(f)
    g=Fun(coefficients(f),S.space)^k
    Fun(coefficients(g),JacobiWeight(k*S.α,k*S.β,space(g)))
end

function .*{JW1<:JacobiWeight,JW2<:JacobiWeight}(f::Fun{JW1},g::Fun{JW2})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β
    m=(Fun(f.coefficients,space(f).space).*Fun(g.coefficients,space(g).space))
    if isapprox(fα+gα,0)&&isapprox(fβ+gβ,0)
        m
    else
        Fun(m.coefficients,JacobiWeight(fα+gα,fβ+gβ,space(m)))
    end
end


./{JW1<:JacobiWeight,JW2<:JacobiWeight}(f::Fun{JW1},g::Fun{JW2})=f*(1/g)

# O(min(m,n)) Ultraspherical conjugated inner product

function conjugatedinnerproduct{λ,D,S,V}(::Type{Ultraspherical{λ,D}},u::Vector{S},v::Vector{V})
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

function conjugatedinnerproduct{C<:Chebyshev}(::Type{C},u::Vector,v::Vector)
    mn = min(length(u),length(v))
    if mn > 1
        return (2u[1]*v[1]+dotu(u[2:mn],v[2:mn]))*π/2
    elseif mn > 0
        return u[1]*v[1]*π
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function conjugatedinnerproduct{D}(::Type{Ultraspherical{1,D}},u::Vector,v::Vector)
    mn = min(length(u),length(v))
    if mn > 0
        return dotu(u[1:mn],v[1:mn])*π/2
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function dotu{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D},D}},g::Fun{Ultraspherical{λ,D}})
    @assert domain(f) == domain(g)
    if f.space.α == f.space.β == λ-0.5
        return complexlength(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdotu(f,g)
    end
end

function dotu{λ,D}(f::Fun{Ultraspherical{λ,D}},g::Fun{JacobiWeight{Ultraspherical{λ,D},D}})
    @assert domain(f) == domain(g)
    if g.space.α == g.space.β == λ-0.5
        return complexlength(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdotu(f,g)
    end
end

function dotu{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D},D}},g::Fun{JacobiWeight{Ultraspherical{λ,D},D}})
    @assert domain(f) == domain(g)
    if f.space.α+g.space.α == f.space.β+g.space.β == λ-0.5
        return complexlength(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdotu(f,g)
    end
end

function linedotu{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D},D}},g::Fun{Ultraspherical{λ,D}})
    @assert domain(f) == domain(g)
    if f.space.α == f.space.β == λ-0.5
        return length(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultlinedotu(f,g)
    end
end

function linedotu{λ,D}(f::Fun{Ultraspherical{λ,D}},g::Fun{JacobiWeight{Ultraspherical{λ,D},D}})
    @assert domain(f) == domain(g)
    if g.space.α == g.space.β == λ-0.5
        return length(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultlinedotu(f,g)
    end
end

function linedotu{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D},D}},g::Fun{JacobiWeight{Ultraspherical{λ,D},D}})
    @assert domain(f) == domain(g)
    if f.space.α+g.space.α == f.space.β+g.space.β == λ-0.5
        return length(domain(f))/2*conjugatedinnerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultlinedotu(f,g)
    end
end
