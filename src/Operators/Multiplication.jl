export Multiplication

abstract AbstractMultiplication{T} <:BandedOperator{T}

immutable Multiplication{D<:FunctionSpace,S<:FunctionSpace,T,V} <: AbstractMultiplication{V}
    f::Fun{D,T}
    space::S

    Multiplication(f::Fun{D,T},sp::S)=new(f,sp)
end

Multiplication{D,T}(f::Fun{D,T},sp::FunctionSpace{ComplexBasis})=Multiplication{D,typeof(sp),T,promote_type(T,Complex{real(T)})}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)
Multiplication{D,T}(f::Fun{D,T},sp::FunctionSpace)=Multiplication{D,typeof(sp),T,T}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)

Multiplication(f::Fun)=Multiplication(f,UnsetSpace())
Multiplication(c::Number)=ConstantOperator(c)

# This covers right multiplication unless otherwise specified.
Multiplication{D,T}(S::FunctionSpace,f::Fun{D,T}) = Multiplication(f,S)


Base.convert{BT<:Operator,S,V,T}(::Type{BT},C::Multiplication{S,V,T})=Multiplication{S,V,T,eltype(BT)}(C.f,C.space)

domainspace{D,S,T,V}(M::Multiplication{D,S,T,V})=M.space
domain(T::Multiplication)=domain(T.f)


## Default implementation: try converting to space of M.f

rangespace{F,T}(D::Multiplication{F,UnsetSpace,T})=UnsetSpace()
bandinds{F,T}(D::Multiplication{F,UnsetSpace,T})=error("No range space attached to Multiplication")
addentries!{F,T}(D::Multiplication{F,UnsetSpace,T},A,kr)=error("No range space attached to Multiplication")


function addentries!{F,S,T}(D::Multiplication{F,S,T},A,kr)
    # Default is to convert to space of f
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override addentries! on Multiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*") for range type"*string(typeof(kr)))
    end
    addentries!(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]),A,kr)
end

function bandinds{F,S,T}(D::Multiplication{F,S,T})
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override bandinds for Multiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*")")
    end
    bandinds(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]))
end

# corresponds to default implementation
function rangespace{F,S,T}(D::Multiplication{F,S,T})
    sp=domainspace(D)
    csp=space(D.f)
    if csp==sp
        error("Override rangespace for Multiplication(::Fun{"*string(typeof(space(D.f)))*",T},"*string(typeof(sp))*")")
    end
    rangespace(TimesOperator([Multiplication(D.f,csp),Conversion(sp,csp)]))
end






##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::AbstractMultiplication,sp::UnsetSpace)=D
promotedomainspace(D::AbstractMultiplication,sp::AnySpace)=D
promotedomainspace(D::AbstractMultiplication,sp::FunctionSpace)=Multiplication(D.f,sp)

choosedomainspace{D}(M::Multiplication{D,UnsetSpace},::AnySpace)=space(M.f)
choosedomainspace{D}(M::Multiplication{D,UnsetSpace},sp)=sp  # we assume multiplication maps spaces to themselves


Base.diagm(a::Fun)=Multiplication(a)

immutable MultiplicationWrapper{D<:FunctionSpace,O<:BandedOperator,V<:Number,T<:Number} <: AbstractMultiplication{T}
    f::Fun{D,V}
    op::O
end

MultiplicationWrapper{D<:FunctionSpace,V<:Number,T<:Number}(f::Fun{D,V},op::BandedOperator{T})=MultiplicationWrapper{D,typeof(op),V,T}(f,op)

addentries!(D::MultiplicationWrapper,A,k::Range)=addentries!(D.op,A,k)
for func in (:rangespace,:domainspace,:bandinds,:domain,:(Base.stride))
    @eval $func(D::MultiplicationWrapper)=$func(D.op)
end

Base.convert{BT<:Operator,S,V,VV,T}(::Type{BT},C::MultiplicationWrapper{S,V,VV,T})=MultiplicationWrapper{S,V,VV,eltype(BT)}(C.f,C.op)




## Multiplication operators allow us to multiply two spaces

# Overrideable
# This should be overriden whenever the multiplication space is different
function .*{T,N,S,V}(f::Fun{S,T},g::Fun{V,N})
    # When the spaces differ we promote and multiply
    if domainscompatible(space(f),space(g))
        # THe bandwidth of Mutliplication is
        # usually the length of the function
        if length(f)≤length(g)
            Multiplication(f,space(g))*g
        else
            Multiplication(g,space(f))*f
        end
    else
        sp=union(space(f),space(g))
        Fun(f,sp).*Fun(g,sp)
    end
end


function transformtimes(f::Fun,g::Fun,n)
    @assert spacescompatible(space(f),space(g))
    f2 = pad(f,n); g2 = pad(g,n)

    sp=space(f)
    chop!(Fun(transform(sp,values(f2).*values(g2)),sp),10eps())
end
transformtimes(f::Fun,g::Fun)=transformtimes(f,g,length(f) + length(g) - 1)
