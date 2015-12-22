export LaurentDirichlet

##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
# TODO: restructure as SumSpace made out of HardyDirichlet
##



immutable LaurentDirichlet{DD} <: UnivariateSpace{ComplexBasis,DD}
    domain::DD
    LaurentDirichlet(d::DD)=new(d)
end
LaurentDirichlet{T<:Number}(d::Vector{T})=LaurentDirichlet(PeriodicDomain(d))
LaurentDirichlet(d::Domain)=LaurentDirichlet{typeof(d)}(d)
LaurentDirichlet()=LaurentDirichlet(PeriodicInterval())

spacescompatible(a::LaurentDirichlet,b::LaurentDirichlet)=domainscompatible(a,b)

canonicalspace(S::LaurentDirichlet)=Laurent(domain(S))

bandinds{DD}(::ConcreteConversion{LaurentDirichlet{DD},Laurent{DD}})=0,2
function addentries!{DD}(C::ConcreteConversion{LaurentDirichlet{DD},Laurent{DD}},A,kr::Range,::Colon)
    if 1 in kr
        A[1,2]+=1
    end
    toeplitz_addentries!([],[1.,0.,1.],A,kr)
end


conversion_rule{DD}(b::LaurentDirichlet,a::Laurent{DD})=b

differentiate{DD}(f::Fun{LaurentDirichlet{DD}}) = differentiate(Fun(f,Laurent))

.*{DD}(f::Fun{Laurent{DD}},g::Fun{LaurentDirichlet{DD}}) = f.*Fun(g,Laurent)
.*{DD}(f::Fun{LaurentDirichlet{DD}},g::Fun{Laurent{DD}}) = Fun(f,Laurent).*g
.*{DD}(f::Fun{LaurentDirichlet{DD}},g::Fun{LaurentDirichlet{DD}}) = Fun(f,Laurent).*g

Base.real{DD}(f::Fun{LaurentDirichlet{DD}}) = real(Fun(f,Laurent))
Base.imag{DD}(f::Fun{LaurentDirichlet{DD}}) = imag(Fun(f,Laurent))


coefficients(v::Vector,::Laurent,::LaurentDirichlet)=laurentdirichlettransform!(copy(v))
coefficients(v::Vector,::LaurentDirichlet,::Laurent)=laurentidirichlettransform!(copy(v))

function laurentdirichlettransform!(w::Vector)
    if length(w) > 1
        for k=length(w)-2:-1:1
            @inbounds w[k] -= w[k+2]
        end
        @inbounds w[1] -= w[2]
    end

    w
end

function laurentidirichlettransform!(w::Vector)
    if length(w) == 1
    elseif length(w) == 2
        @inbounds w[1] += w[2]
    else
        @inbounds w[1] += w[2]+w[3]
        for k=4:length(w)
            @inbounds w[k-2] += w[k]
        end
    end

    w
end


##
# CosDirichlet represents
# 1,cos(θ)+1,cos(2θ)+cos(θ),…
#


immutable CosDirichlet{DD} <: RealUnivariateSpace{DD}
    domain::DD
end

CosDirichlet()=CosDirichlet(PeriodicInterval())

spacescompatible(a::CosDirichlet,b::CosDirichlet)=domainscompatible(a,b)

canonicalspace(S::CosDirichlet)=CosSpace(domain(S))

bandinds{CS<:CosSpace,CD<:CosDirichlet}(::ConcreteConversion{CD,CS})=0,1
addentries!{CS<:CosSpace,CD<:CosDirichlet}(C::ConcreteConversion{CD,CS},A,kr::Range,::Colon)=toeplitz_addentries!([],[1.,1.],A,kr)


conversion_rule(b::CosDirichlet,a::CosSpace)=b
