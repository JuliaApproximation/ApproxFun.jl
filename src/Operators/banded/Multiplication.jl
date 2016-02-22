export Multiplication

abstract Multiplication{T} <:BandedOperator{T}

immutable ConcreteMultiplication{D<:Space,S<:Space,V,T} <: Multiplication{T}
    f::Fun{D,V}
    space::S

    ConcreteMultiplication(f::Fun{D,V},sp::S)=new(f,sp)
end

function ConcreteMultiplication{D,T}(f::Fun{D,T},sp::Space)
    @assert domainscompatible(space(f),sp)
    ConcreteMultiplication{D,typeof(sp),T,mat_promote_type(T,eltype(sp))}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)
end


Multiplication(f::Fun,sp::Space)=ConcreteMultiplication(f,sp)
Multiplication(f::Fun)=Multiplication(f,UnsetSpace())
Multiplication(c::Number)=Multiplication(Fun(c) )

# This covers right multiplication unless otherwise specified.
Multiplication{D,T}(S::Space,f::Fun{D,T}) = Multiplication(f,S)

for TYP in (:Operator,:BandedOperator)
    @eval function Base.convert{S,V,TT,T}(::Type{$TYP{T}},C::ConcreteMultiplication{S,V,TT})
        if T==eltype(C)
            C
        else
            ConcreteMultiplication{S,V,TT,T}(C.f,C.space)
        end
    end
end

domainspace{D,S,T,V}(M::ConcreteMultiplication{D,S,T,V})=M.space
domain(T::ConcreteMultiplication)=domain(T.f)


## Default implementation: try converting to space of M.f

# avoid ambiguity
rangespace{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T})=UnsetSpace()
bandinds{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T})=error("No range space attached to Multiplication")
addentries!{F,T}(D::ConcreteMultiplication{F,UnsetSpace,T},A,kr,::Colon)=error("No range space attached to Multiplication")


function addentries!{F,S,T}(D::ConcreteMultiplication{F,S,T},A,kr,::Colon)
    # Default is to convert to space of f
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override addentries! on ConcreteMultiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*") for range type"*string(typeof(kr)))
    end
    addentries!(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]),A,kr,:)
end

function bandinds{F,S,T}(D::ConcreteMultiplication{F,S,T})
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override bandinds for ConcreteMultiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*")")
    end
    bandinds(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]))
end

# corresponds to default implementation
function rangespace{F,S,T}(D::ConcreteMultiplication{F,S,T})
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override rangespace for ConcreteMultiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*")")
    end
    rangespace(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]))
end






##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::Multiplication,sp::UnsetSpace)=D
promotedomainspace(D::Multiplication,sp::AnySpace)=D
promotedomainspace(D::Multiplication,sp::Space)=Multiplication(D.f,sp)

choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},::AnySpace)=space(M.f)
choosedomainspace{D}(M::ConcreteMultiplication{D,UnsetSpace},sp)=sp  # we assume multiplication maps spaces to themselves


Base.diagm(a::Fun)=Multiplication(a)

immutable MultiplicationWrapper{D<:Space,O<:BandedOperator,V,T} <: Multiplication{T}
    f::Fun{D,V}
    op::O
end

MultiplicationWrapper{D<:Space,V}(T::Type,f::Fun{D,V},op::BandedOperator)=MultiplicationWrapper{D,typeof(op),V,T}(f,op)
MultiplicationWrapper{D<:Space,V}(f::Fun{D,V},op::BandedOperator)=MultiplicationWrapper(eltype(op),f,op)

addentries!(D::MultiplicationWrapper,A,k::Range,::Colon)=addentries!(D.op,A,k,:)
for func in (:rangespace,:domainspace,:bandinds,:domain,:(Base.stride))
    @eval $func(D::MultiplicationWrapper)=$func(D.op)
end

Base.convert{BT<:MultiplicationWrapper}(::Type{BT},C::BT)=C
function Base.convert{BT<:Operator,S,V,VV,T}(::Type{BT},C::MultiplicationWrapper{S,V,VV,T})
    if eltype(BT)==eltype(C)
        C
    else
        MultiplicationWrapper{S,V,VV,eltype(BT)}(C.f,C.op)
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
        m,n = length(f),length(g)
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
transformtimes(f::Fun,g::Fun)=transformtimes(f,g,length(f) + length(g) - 1)
