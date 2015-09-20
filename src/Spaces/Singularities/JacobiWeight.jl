


export JacobiWeight


##
# WeightSpace represents a space that weights another space
##

abstract WeightSpace <: RealUnivariateSpace{Interval{Float64}}  #TODO: Why Interval?


domain(S::WeightSpace)=domain(S.space)


points(sp::WeightSpace,n)=points(sp.space,n)
plan_transform(S::WeightSpace,vals::Vector)=1./weight(S,points(S,length(vals))),plan_transform(S.space,vals)
plan_itransform(S::WeightSpace,vals::Vector)=weight(S,points(S,length(vals))),plan_itransform(S.space,vals)


transform(sp::WeightSpace,vals::Vector,plan)=transform(sp.space,vals.*plan[1],plan[2])
itransform(sp::WeightSpace,cfs::Vector,plan)=itransform(sp.space,cfs,plan[2]).*plan[1]



function evaluate{WS<:WeightSpace,T}(f::Fun{WS,T},x)
    tol=1.0E-14
    fv=Fun(f.coefficients,space(f).space)[x]
    if isa(fv,Number)&&abs(fv)<tol
        #TODO: Why this special case??
        zero(T)
    else
        weight(space(f),x).*fv
    end
end


##
# JacobiWeight
# represents a function on [-1,1] weighted by (1+x)^α*(1-x)^β
# note the inconsistency of the parameters with Jacobi
# when the domain is [a,b] the weight is inferred by mapping to [-1,1]
##


immutable JacobiWeight{S} <: WeightSpace
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

JacobiWeight(a::Number,b::Number,d::Space)=JacobiWeight{typeof(d)}(@compat(Float64(a)),@compat(Float64(b)),d)
JacobiWeight(a::Number,b::Number,d::IntervalDomain)=JacobiWeight(@compat(Float64(a)),@compat(Float64(b)),Space(d))
JacobiWeight(a::Number,b::Number,d::Vector)=JacobiWeight(@compat(Float64(a)),@compat(Float64(b)),Space(d))
JacobiWeight(a::Number,b::Number)=JacobiWeight(a,b,Chebyshev())

JacobiWeight(a::Number,b::Number,s::Vector) = map(s->JacobiWeight(a,b,s),s)
JacobiWeight(a::Number,b::Number,s::PiecewiseSpace) = PiecewiseSpace(JacobiWeight(a,b,vec(s)))

identity_fun(S::JacobiWeight)=isapproxinteger(S.α)&&isapproxinteger(S.β)?Fun(x->x,S):Fun(identity,domain(S))



spacescompatible(A::JacobiWeight,B::JacobiWeight)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)
spacescompatible{T,D<:Interval}(A::JacobiWeight,B::Space{T,D})=spacescompatible(A,JacobiWeight(0,0,B))
spacescompatible{T,D<:Interval}(B::Space{T,D},A::JacobiWeight)=spacescompatible(A,JacobiWeight(0,0,B))

transformtimes{S,V}(f::Fun{JacobiWeight{S}},g::Fun{JacobiWeight{V}}) = Fun(coefficients(transformtimes(Fun(f.coefficients,f.space.space),Fun(g.coefficients,g.space.space))),JacobiWeight(f.space.α+g.space.α,f.space.β+g.space.β,f.space.space))
transformtimes{S}(f::Fun{JacobiWeight{S}},g::Fun) = Fun(coefficients(transformtimes(Fun(f.coefficients,f.space.space),g)),f.space)
transformtimes{S}(f::Fun,g::Fun{JacobiWeight{S}}) = Fun(coefficients(transformtimes(Fun(g.coefficients,g.space.space),f)),g.space)

##  α and β are opposite the convention for Jacobi polynomials
# Here, α is the left algebraic singularity and β is the right algebraic singularity.

jacobiweight(α,β,x)=(1+x).^α.*(1-x).^β
weight(sp::JacobiWeight,x)=jacobiweight(sp.α,sp.β,tocanonical(sp,x))



setdomain(sp::JacobiWeight,d::Domain)=JacobiWeight(sp.α,sp.β,setdomain(sp.space,d))

# we assume that points avoids singularities


##TODO: paradigm for same space
function coefficients(f::Vector,sp1::JacobiWeight,sp2::JacobiWeight)
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β

    if isapprox(c,α) && isapprox(d,β)
        coefficients(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end
coefficients{S,n,st}(f::Vector,sp::JacobiWeight,S2::SliceSpace{n,st,S,RealBasis})=error("Implement")
coefficients{S,n,st}(f::Vector,S2::SliceSpace{n,st,S,RealBasis},sp::JacobiWeight)=error("Implement")
#TODO: it could be possible that we want to JacobiWeight a SumSpace....
coefficients{SV,T,D<:Interval}(f::Vector,sp::JacobiWeight,S2::SumSpace{SV,T,D,1})=sumspacecoefficients(f,sp,S2)
coefficients{T,D<:Interval}(f::Vector,sp::JacobiWeight,S2::Space{T,D,1})=coefficients(f,sp,JacobiWeight(0,0,S2))
coefficients{T,D<:Interval}(f::Vector,S2::Space{T,D,1},sp::JacobiWeight)=coefficients(f,JacobiWeight(0,0,S2),sp)

increase_jacobi_parameter(f)=Fun(f,JacobiWeight(f.space.α+1,f.space.β+1,space(f).space))
increase_jacobi_parameter(s,f)=s==-1?Fun(f,JacobiWeight(f.space.α+1,f.space.β,space(f).space)):Fun(f,JacobiWeight(f.space.α,f.space.β+1,space(f).space))



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
        function ($op){S}(c::Number,f::Fun{JacobiWeight{S}})
            g=($op)(c,Fun(f.coefficients,space(f).space))
            Fun(g.coefficients,JacobiWeight(-f.space.α,-f.space.β,space(g)))
        end
    end
end

function .^{J<:JacobiWeight}(f::Fun{J},k::Float64)
    S=space(f)
    g=Fun(coefficients(f),S.space)^k
    Fun(coefficients(g),JacobiWeight(k*S.α,k*S.β,space(g)))
end

function .*{S,V}(f::Fun{JacobiWeight{S}},g::Fun{JacobiWeight{V}})
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


./{T,N}(f::Fun{JacobiWeight{T}},g::Fun{JacobiWeight{N}})=f*(1/g)

# O(min(m,n)) Ultraspherical inner product

function innerproduct{λ,D,S,V}(::Type{Ultraspherical{λ,D}},u::Vector{S},v::Vector{V})
    T,mn = promote_type(S,V),min(length(u),length(v))
    if mn > 1
        wi = sqrt(convert(T,π))*gamma(λ+one(T)/2)/gamma(λ+one(T))
        ret = conj(u[1])*wi*v[1]
        for i=2:mn
          wi *= (i-2one(T)+2λ)/(i-one(T)+λ)*(i-2one(T)+λ)/(i-one(T))
          ret += conj(u[i])*wi*v[i]
        end
        return ret
    elseif mn > 0
        wi = sqrt(convert(T,π))*gamma(λ+one(T)/2)/gamma(λ+one(T))
        return conj(u[1])*wi*v[1]
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function innerproduct{C<:Chebyshev}(::Type{C},u::Vector,v::Vector)
    mn = min(length(u),length(v))
    if mn > 1
        return (2conj(u[1])*v[1]+dot(u[2:mn],v[2:mn]))*π/2
    elseif mn > 0
        return conj(u[1])*v[1]*π
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function innerproduct{D}(::Type{Ultraspherical{1,D}},u::Vector,v::Vector)
    mn = min(length(u),length(v))
    if m > 0
        return dot(u[1:mn],v[1:mn])*π/2
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function Base.dot{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D}}},g::Fun{Ultraspherical{λ,D}})
    @assert domain(f) == domain(g)
    if f.space.α == f.space.β == λ-0.5
        return complexlength(domain(f))/2*innerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdot(f,g)
    end
end

function Base.dot{λ,D}(f::Fun{Ultraspherical{λ,D}},g::Fun{JacobiWeight{Ultraspherical{λ,D}}})
    @assert domain(f) == domain(g)
    if g.space.α == g.space.β == λ-0.5
        return complexlength(domain(f))/2*innerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdot(f,g)
    end
end

function Base.dot{λ,D}(f::Fun{JacobiWeight{Ultraspherical{λ,D}}},g::Fun{JacobiWeight{Ultraspherical{λ,D}}})
    @assert domain(f) == domain(g)
    if f.space.α+g.space.α == f.space.β+g.space.β == λ-0.5
        return complexlength(domain(f))/2*innerproduct(Ultraspherical{λ,D},f.coefficients,g.coefficients)
    else
        return defaultdot(f,g)
    end
end
