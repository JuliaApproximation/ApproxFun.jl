type ProductDomain{D<:Domain}
    domains::Vector{D} 
end

ProductDomain(A,B)=ProductDomain([A,B])
*(A::Domain,B::Domain)=ProductDomain(A,B)

domain(f::Fun2D)=domain(f.A[1])*domain(f.B[1])

Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]