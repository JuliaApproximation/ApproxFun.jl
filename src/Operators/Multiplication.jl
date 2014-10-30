abstract AbstractMultiplication{T} <:BandedOperator{T}

immutable Multiplication{D<:FunctionSpace,S<:FunctionSpace,T<:Number} <: AbstractMultiplication{T}
    f::Fun{D,T}
    space::S
end

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


immutable MultiplicationWrapper{O<:BandedOperator,T<:Number} <: AbstractMultiplication{T}
    op::O
end

MultiplicationWrapper{T<:Number}(op::BandedOperator{T})=MultiplicationWrapper{typeof(op),T}(op)

addentries!(D::MultiplicationWrapper,A::ShiftArray,k::Range)=addentries!(D.op,A,k)
for func in (:rangespace,:domainspace,:bandinds,:domain)
    @eval $func(D::MultiplicationWrapper)=$func(D.op)
end


