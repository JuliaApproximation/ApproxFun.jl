export LaurentDirichlet

##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
# TODO: restructure as SumSpace made out of HardyDirichlet
##



immutable LaurentDirichlet{DD} <: UnivariateSpace{ComplexBasis,DD}
    domain::DD
    (::Type{LaurentDirichlet{DD}}){DD}(d::DD) = new{DD}(d)
end
LaurentDirichlet{T<:Number}(d::Vector{T})=LaurentDirichlet(PeriodicDomain(d))
LaurentDirichlet(d::Domain)=LaurentDirichlet{typeof(d)}(d)
LaurentDirichlet()=LaurentDirichlet(PeriodicInterval())

spacescompatible(a::LaurentDirichlet,b::LaurentDirichlet)=domainscompatible(a,b)

canonicalspace(S::LaurentDirichlet)=Laurent(domain(S))


Conversion{DD}(a::LaurentDirichlet{DD},b::Laurent{DD})=ConcreteConversion(a,b)
bandinds{DD}(::ConcreteConversion{LaurentDirichlet{DD},Laurent{DD}})=0,2

function getindex{DD}(C::ConcreteConversion{LaurentDirichlet{DD},Laurent{DD}},k::Integer,j::Integer)
    if k==1 && j==2
        one(eltype(C))
    elseif k==j || j==k+2
        one(eltype(C))
    else
        zero(eltype(C))
    end
end

conversion_rule{DD}(b::LaurentDirichlet,a::Laurent{DD})=b

differentiate{DD}(f::Fun{LaurentDirichlet{DD}}) = differentiate(Fun(f,Laurent))

for op in (:+,:-,:*)
    @eval begin
        $op{DD}(f::Fun{Laurent{DD}},g::Fun{LaurentDirichlet{DD}}) = $op(f,Fun(g,Laurent))
        $op{DD}(f::Fun{LaurentDirichlet{DD}},g::Fun{Laurent{DD}}) = $op(Fun(f,Laurent),g)
        $op{DD}(f::Fun{LaurentDirichlet{DD}},g::Fun{LaurentDirichlet{DD}}) = $op(Fun(f,Laurent),Fun(g,Laurent))
    end
end


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

spacescompatible(a::CosDirichlet,b::CosDirichlet) = domainscompatible(a,b)

canonicalspace(S::CosDirichlet)=CosSpace(domain(S))

Conversion(a::CosSpace,b::CosDirichlet) = ConcreteConversion(a,b)
bandinds{CS<:CosSpace,CD<:CosDirichlet}(::ConcreteConversion{CD,CS})=0,1
getindex{CS<:CosSpace,CD<:CosDirichlet}(C::ConcreteConversion{CD,CS},k::Integer,j::Integer) =
    (k==j||j==k+1)?one(eltype(C)):zero(eltype(C))


conversion_rule(b::CosDirichlet,a::CosSpace)=b
