abstract AbstractMultiplication{T} <:BandedOperator{T}

immutable Multiplication{D<:FunctionSpace,S<:FunctionSpace,T<:Number} <: AbstractMultiplication{T}
    f::Fun{D,T}
    space::S
    
    Multiplication(f::Fun{D,T},sp::S)=new(f,sp)
end

Multiplication{D,T,S}(f::Fun{D,T},sp::S)=Multiplication{D,S,T}(f,sp)

Multiplication(f::Fun)=Multiplication(f,AnySpace())

Multiplication(c::Number)=ConstantOperator(c)




domainspace{D,S,T}(M::Multiplication{D,S,T})=M.space
rangespace{D,S,T}(M::Multiplication{D,S,T})=M.space





bandinds(T::Multiplication)=(1-length(T.f.coefficients),length(T.f.coefficients)-1)
domain(T::Multiplication)=domain(T.f)



##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::AbstractMultiplication,sp::AnySpace)=D
promotedomainspace(D::AbstractMultiplication,sp::FunctionSpace)=Multiplication(D.f,sp)


Base.diagm(a::Fun)=Multiplication(a)


immutable MultiplicationWrapper{D<:FunctionSpace,O<:BandedOperator,T<:Number} <: AbstractMultiplication{T}
    f::Fun{D,T}
    op::O
end

MultiplicationWrapper{D<:FunctionSpace,T<:Number}(f::Fun{D,T},op::BandedOperator{T})=MultiplicationWrapper{D,typeof(op),T}(f,op)

addentries!(D::MultiplicationWrapper,A,k::Range)=addentries!(D.op,A,k)
for func in (:rangespace,:domainspace,:bandinds,:domain)
    @eval $func(D::MultiplicationWrapper)=$func(D.op)
end





## Multiplication operators allowus to multiply two spaces

# Overrideable
# This should be overriden whenever the multiplication space is different
function .*{T,N,S,V}(f::Fun{S,T},g::Fun{V,N})
    # When the spaces differ we promote and multiply   
    if domainscompatible(space(f),space(g))
        Multiplication(f,space(g))*g
    else         
        sp=union(space(f),space(g))
        Fun(f,sp).*Fun(g,sp)
    end
end