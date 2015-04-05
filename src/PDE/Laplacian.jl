

@calculus_operator(Laplacian,AbstractLaplacian,LaplacianWrapper)

Laplacian(S::BivariateSpace,k)=Laplacian{typeof(S),Int,BandedMatrix{eltype(S)}}(S,k)


immutable Dirichlet{S,T} <: Operator{T}
    space::S
end
Dirichlet(sp::BivariateSpace)=Dirichlet{typeof(sp),BandedMatrix{eltype(sp)}}(sp)
Dirichlet(d::Domain)=Dirichlet(Space(d))

domainspace(S::Dirichlet)=S.space
rangespace(B::Dirichlet)=Space(∂(domain(B)))


