##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
# TODO: restructure as SumSpace made out of HardyDirichlet
##



immutable LaurentDirichlet <: UnivariateSpace{ComplexBasis}
    domain::PeriodicDomain
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


##
# CosDirichlet represents
# 1,cos(θ)+1,cos(2θ)+cos(θ),…
#


immutable CosDirichlet <: RealUnivariateSpace
    domain::PeriodicDomain
end

CosDirichlet()=CosDirichlet(PeriodicInterval())

spacescompatible(a::CosDirichlet,b::CosDirichlet)=domainscompatible(a,b)

canonicalspace(S::CosDirichlet)=CosSpace(domain(S))

bandinds(::Conversion{CosDirichlet,CosSpace})=0,1
addentries!(C::Conversion{CosDirichlet,CosSpace},A,kr::Range)=toeplitz_addentries!([],[1.,1.],A,kr)


conversion_rule(b::CosDirichlet,a::CosSpace)=b

