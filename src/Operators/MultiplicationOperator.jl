export MultiplicationOperator



type MultiplicationOperator{T<:Number,D<:FunctionSpace} <: BandedOperator{T}
    f::Fun{T,D}
end

MultiplicationOperator(c::Number)=ConstantOperator(c)




domainspace(M::MultiplicationOperator)=space(M.f)
rangespace(M::MultiplicationOperator)=space(M.f)





bandinds(T::MultiplicationOperator)=(1-length(T.f.coefficients),length(T.f.coefficients)-1)
domain(T::MultiplicationOperator)=domain(T.f)



##multiplication can always be promoted, range space is allowed to change
promotedomainspace(D::MultiplicationOperator,sp::AnySpace)=D
promotedomainspace(D::MultiplicationOperator,sp::FunctionSpace)=MultiplicationOperator(Fun(D.f,sp))


