export Multiplication

abstract type Multiplication{D,S,T} <:Operator{T} end

struct ConcreteMultiplication{D<:Space,S<:Space,T} <: Multiplication{D,S,T}
    f::VFun{D,T}
    space::S

    ConcreteMultiplication{D,S,T}(f::Fun{D,T},sp::S) where {D,S,T} = new{D,S,T}(f,sp)
end

function ConcreteMultiplication{V,D,T}(::Type{V},f::Fun{D,T},sp::Space)
    if !domainscompatible(space(f),sp)
        error("Domain mismatch: cannot multiply function on $(domain(f)) to function on $(domain(sp))")
    end
    ConcreteMultiplication{D,typeof(sp),V}(
        convert(Fun{D,V},chop(f,maximum(abs,f.coefficients)*40*eps(eltype(f)))),sp)
end


function ConcreteMultiplication{D,T}(f::Fun{D,T},sp::Space)
    if !domainscompatible(space(f),sp)
        error("Domain mismatch: cannot multiply function on $(domain(f)) to function on $(domain(sp))")
    end
    V = promote_type(T,eltype(sp))
    ConcreteMultiplication{D,typeof(sp),V}(convert(Fun{D,V},chop(f,maximum(abs,f.coefficients)*40*eps(eltype(f)))),sp)
end

# We do this in two stages to support Modifier spaces
# without ambiguity errors
function defaultMultiplication(f::Fun,sp::Space)
    csp=space(f)
    if csp==sp
        error("Implement Multiplication(::Fun{$(typeof(space(f)))},::$(typeof(sp)))")
    end
    MultiplicationWrapper(f,Multiplication(f,csp)*Conversion(sp,csp))
end

Multiplication(f::Fun,sp::Space) = defaultMultiplication(f,sp)


Multiplication(f::Fun,sp::UnsetSpace) = ConcreteMultiplication(f,sp)
Multiplication(f::Fun) = Multiplication(f,UnsetSpace())

Multiplication(c::Number,sp::Space) = Multiplication(Fun(c),sp)
Multiplication(sp::Space,c::Number) = Multiplication(sp,Fun(c))
Multiplication(c::Number) = Multiplication(Fun(c) )

# This covers right multiplication unless otherwise specified.
Multiplication(S::Space,f::Fun) = Multiplication(f,S)


function Base.convert{S,V,T}(::Type{Operator{T}},C::ConcreteMultiplication{S,V})
    if T==eltype(C)
        C
    else
        ConcreteMultiplication{S,V,T}(Fun{S,T}(C.f),C.space)
    end
end

domainspace{D,S,T}(M::ConcreteMultiplication{D,S,T}) = M.space
domain(T::ConcreteMultiplication) = domain(T.f)


## Default implementation: try converting to space of M.f

# avoid ambiguity
rangespace{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T}) = UnsetSpace()
getindex{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T},k::Integer,j::Integer) =
    error("No range space attached to Multiplication")






##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::Multiplication,sp::UnsetSpace) = D
promotedomainspace(D::Multiplication,sp::Space) = Multiplication(D.f,sp)
promoterangespace{P}(D::ConcreteMultiplication{P,UnsetSpace},sp::UnsetSpace) = D
promoterangespace{P}(D::ConcreteMultiplication{P,UnsetSpace},sp::Space) =
    promoterangespace(Multiplication(D.f,ConstantSpace(domain(sp))), sp)

choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},::UnsetSpace) = space(M.f)
# we assume multiplication maps spaces to themselves
choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},sp::Space) = sp


Base.diagm(a::Fun) = Multiplication(a)

struct MultiplicationWrapper{D<:Space,S<:Space,O<:Operator,T} <: Multiplication{D,S,T}
    f::VFun{D,T}
    op::O
end

MultiplicationWrapper{D<:Space,V}(T::Type,f::Fun{D,V},op::Operator) = MultiplicationWrapper{D,typeof(domainspace(op)),typeof(op),T}(f,op)
MultiplicationWrapper{D<:Space,V}(f::Fun{D,V},op::Operator) = MultiplicationWrapper(eltype(op),f,op)

@wrapper MultiplicationWrapper

function Base.convert{TT,S,V,O,T}(::Type{Operator{TT}},C::MultiplicationWrapper{S,V,O,T})
    if TT==T
        C
    else
        MultiplicationWrapper(Fun{S,TT}(C.f),Operator{TT}(C.op))::Operator{TT}
    end
end


## Multiplication operators allow us to multiply two spaces



hasfasttransform(::)=false
hasfasttransform(f::Fun)=hasfasttransform(space(f))
hasfasttransformtimes(f,g)=spacescompatible(f,g) && hasfasttransform(f) && hasfasttransform(g)


# This should be overriden whenever the multiplication space is different
function *{T,N,S,V}(f::Fun{S,T},g::Fun{V,N})
    # When the spaces differ we promote and multiply
    if domainscompatible(space(f),space(g))
        m,n = ncoefficients(f),ncoefficients(g)
        # Heuristic division of parameter space between value-space and coefficient-space multiplication.
        if hasfasttransformtimes(f,g) && log10(m)*log10(n)>4
            transformtimes(f,g)
        elseif m≤n
            coefficienttimes(f,g)
        else
            coefficienttimes(g,f)
        end
    else
        sp=union(space(f),space(g))
        Fun(f,sp)*Fun(g,sp)
    end
end

coefficienttimes(f::Fun,g::Fun) = Multiplication(f,space(g))*g

function transformtimes(f::Fun,g::Fun,n)
    @assert spacescompatible(space(f),space(g))
    isempty(f.coefficients) && return f
    isempty(g.coefficients) && return g
    f2,g2,sp = pad(f,n),pad(g,n),space(f)
    hc = transform(sp,values(f2).*values(g2))
    chop!(Fun(sp,hc),10eps(eltype(hc)))
end
transformtimes(f::Fun,g::Fun) = transformtimes(f,g,ncoefficients(f) + ncoefficients(g) - 1)



*(a::Fun,L::UniformScaling) = Multiplication(a*L.λ,UnsetSpace())
*(L::UniformScaling,a::Fun) = L.λ*a


## docs

doc"""
`Multiplication(f::Fun,sp::Space)` is the operator representing multiplication by
`f` on functions in the space `sp`.
"""
Multiplication(::Fun,::Space)

doc"""
`Multiplication(f::Fun)` is the operator representing multiplication by
`f` on an unset space of functions.  Spaces will be inferred when applying or
manipulating the operator.
"""
Multiplication(::Fun)
