immutable ProductDomain{D<:Domain}
    domains::Vector{D} 
end

ProductDomain(A,B)=ProductDomain([A,B])
*(A::Domain,B::Domain)=ProductDomain(A,B)


Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]


immutable TensorSpace{S<:FunctionSpace,T<:FunctionSpace}
    spaces::(S,T)
end

TensorSpace(A,B)=TensorSpace((A,B))
⊗(A::FunctionSpace,B::FunctionSpace)=TensorSpace(A,B)


Base.getindex(d::TensorSpace,k::Integer)=d.spaces[k]
