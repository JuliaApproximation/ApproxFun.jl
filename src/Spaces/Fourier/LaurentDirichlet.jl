##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
##


immutable LaurentDirichlet <: FunctionSpace{ComplexBasis}
    domain::Union(PeriodicDomain,AnyDomain)
end

LaurentDirichlet()=LaurentDirichlet(PeriodicInterval())

spacescompatible(a::LaurentDirichlet,b::LaurentDirichlet)=domainscompatible(a,b)

canonicalspace(S::LaurentDirichlet)=Laurent(domain(S))

bandinds(::Conversion{LaurentDirichlet,Laurent})=0,2
function addentries!(C::Conversion{LaurentDirichlet,Laurent},A,kr::Range)
    A[1,2]+=1
    toeplitz_addentries!([],[1.,0.,1.],A,kr)
end


conversion_rule(b::LaurentDirichlet,a::Laurent)=b

