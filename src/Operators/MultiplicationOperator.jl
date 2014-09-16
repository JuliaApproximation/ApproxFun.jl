export MultiplicationOperator



type MultiplicationOperator{T<:Number,D<:FunctionSpace,S<:FunctionSpace} <: BandedOperator{T}
    f::Fun{T,D}
    space::S
end

MultiplicationOperator(f::Fun)=MultiplicationOperator(f,space(f))

MultiplicationOperator(c::Number)=ConstantOperator(c)




domainspace(M::MultiplicationOperator)=M.space
rangespace(M::MultiplicationOperator)=M.space





bandinds(T::MultiplicationOperator)=(1-length(T.f.coefficients),length(T.f.coefficients)-1)
domain(T::MultiplicationOperator)=domain(T.f)



##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::MultiplicationOperator,sp::AnySpace)=D
promotedomainspace(D::MultiplicationOperator,sp::FunctionSpace)=MultiplicationOperator(D.f,sp)


