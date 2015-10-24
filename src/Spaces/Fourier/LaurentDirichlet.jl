##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
# TODO: restructure as SumSpace made out of HardyDirichlet
##



immutable LaurentDirichlet <: UnivariateSpace{ComplexBasis,AnyDomain}
    domain::PeriodicDomain
end

LaurentDirichlet()=LaurentDirichlet(PeriodicInterval())

spacescompatible(a::LaurentDirichlet,b::LaurentDirichlet)=domainscompatible(a,b)

canonicalspace(S::LaurentDirichlet)=Laurent(domain(S))

bandinds{DD}(::Conversion{LaurentDirichlet,Laurent{DD}})=0,2
function addentries!{DD}(C::Conversion{LaurentDirichlet,Laurent{DD}},A,kr::Range,::Colon)
    A[1,2]+=1
    toeplitz_addentries!([],[1.,0.,1.],A,kr)
end


conversion_rule{DD}(b::LaurentDirichlet,a::Laurent{DD})=b


##
# CosDirichlet represents
# 1,cos(θ)+1,cos(2θ)+cos(θ),…
#


immutable CosDirichlet <: RealUnivariateSpace{AnyDomain}
    domain::PeriodicDomain
end

CosDirichlet()=CosDirichlet(PeriodicInterval())

spacescompatible(a::CosDirichlet,b::CosDirichlet)=domainscompatible(a,b)

canonicalspace(S::CosDirichlet)=CosSpace(domain(S))

bandinds{CS<:CosSpace}(::Conversion{CosDirichlet,CS})=0,1
addentries!{CS<:CosSpace}(C::Conversion{CosDirichlet,CS},A,kr::Range,::Colon)=toeplitz_addentries!([],[1.,1.],A,kr)


conversion_rule(b::CosDirichlet,a::CosSpace)=b
