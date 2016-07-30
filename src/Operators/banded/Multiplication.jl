export Multiplication

abstract Multiplication{T} <:Operator{T}

immutable ConcreteMultiplication{D<:Space,S<:Space,V,T} <: Multiplication{T}
    f::Fun{D,V}
    space::S

    ConcreteMultiplication(f::Fun{D,V},sp::S)=new(f,sp)
end

function ConcreteMultiplication{D,T}(f::Fun{D,T},sp::Space)
    @assert domainscompatible(space(f),sp)
    ConcreteMultiplication{D,typeof(sp),T,mat_promote_type(T,eltype(sp))}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)
end

#TODO: SPECIALOPS Remove
isbanded(A::ConcreteMultiplication) = true

# We do this in two stages to support Modifier spaces
# without ambiguity errors
function defaultMultiplication(f::Fun,sp::Space)
    csp=space(f)
    if csp==sp
        error("Implement Multiplication(::Fun{$(typeof(space(f)))},::$(typeof(sp)))")
    end
    MultiplicationWrapper(f,Multiplication(f,csp)*Conversion(sp,csp))
end

Multiplication(f::Fun,sp::Space)=defaultMultiplication(f,sp)


Multiplication(f::Fun,sp::UnsetSpace)=ConcreteMultiplication(f,sp)
Multiplication(f::Fun)=Multiplication(f,UnsetSpace())
Multiplication(c::Number)=Multiplication(Fun(c) )

# This covers right multiplication unless otherwise specified.
Multiplication{D,T}(S::Space,f::Fun{D,T}) = Multiplication(f,S)


function Base.convert{S,V,TT,T}(::Type{Operator{T}},C::ConcreteMultiplication{S,V,TT})
    if T==eltype(C)
        C
    else
        ConcreteMultiplication{S,V,TT,T}(C.f,C.space)
    end
end

domainspace{D,S,T,V}(M::ConcreteMultiplication{D,S,T,V})=M.space
domain(T::ConcreteMultiplication)=domain(T.f)


## Default implementation: try converting to space of M.f

# avoid ambiguity
rangespace{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T}) = UnsetSpace()
getindex{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T},k::Integer,j::Integer) =
    error("No range space attached to Multiplication")






##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::Multiplication,sp::UnsetSpace) = D
promotedomainspace(D::Multiplication,sp::Space) = Multiplication(D.f,sp)

choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},::UnsetSpace) = space(M.f)
# we assume multiplication maps spaces to themselves
choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},sp) = sp


Base.diagm(a::Fun) = Multiplication(a)

immutable MultiplicationWrapper{D<:Space,O<:Operator,V,T} <: Multiplication{T}
    f::Fun{D,V}
    op::O
end

MultiplicationWrapper{D<:Space,V}(T::Type,f::Fun{D,V},op::Operator)=MultiplicationWrapper{D,typeof(op),V,T}(f,op)
MultiplicationWrapper{D<:Space,V}(f::Fun{D,V},op::Operator)=MultiplicationWrapper(eltype(op),f,op)

@wrapper MultiplicationWrapper

function Base.convert{TT,S,V,VV,T}(::Type{Operator{TT}},C::MultiplicationWrapper{S,V,VV,T})
    if TT==eltype(C)
        C
    else
        MultiplicationWrapper{S,V,VV,TT}(C.f,C.op)
    end
end


## Multiplication operators allow us to multiply two spaces



hasfasttransform(::)=false
hasfasttransform(f::Fun)=hasfasttransform(space(f))
hasfasttransformtimes(f,g)=spacescompatible(f,g) && hasfasttransform(f) && hasfasttransform(g)


# This should be overriden whenever the multiplication space is different
function .*{T,N,S,V}(f::Fun{S,T},g::Fun{V,N})
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
        Fun(f,sp).*Fun(g,sp)
    end
end

coefficienttimes(f::Fun,g::Fun) = Multiplication(f,space(g))*g

function transformtimes(f::Fun,g::Fun,n)
    @assert spacescompatible(space(f),space(g))
    f2,g2,sp = pad(f,n),pad(g,n),space(f)
    hc = transform(sp,values(f2).*values(g2))
    chop!(Fun(hc,sp),10norm(hc,Inf)*eps(eltype(hc)))
end
transformtimes(f::Fun,g::Fun)=transformtimes(f,g,ncoefficients(f) + ncoefficients(g) - 1)
