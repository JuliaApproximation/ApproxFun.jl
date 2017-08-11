export LaurentDirichlet

##
# LaurentDirichlet represents functions in the basis
#       1, z^(-1)+1,z+1,z^(-2)+z^(-1),z^2+z,…
# which incorporates decay at z = -1
# TODO: restructure as SumSpace made out of HardyDirichlet
##



struct LaurentDirichlet{DD,RR} <: Space{DD,RR}
    domain::DD
    LaurentDirichlet{DD,RR}(d::DD) where {DD,RR} = new(d)
end
LaurentDirichlet(d::Domain) = LaurentDirichlet{typeof(d),complex(eltype(d))}(d)
LaurentDirichlet() = LaurentDirichlet(PeriodicInterval())

spacescompatible(a::LaurentDirichlet,b::LaurentDirichlet) = domainscompatible(a,b)

canonicalspace(S::LaurentDirichlet) = Laurent(domain(S))


Conversion{DD,RR}(a::LaurentDirichlet{DD,RR},b::Laurent{DD,RR}) = ConcreteConversion(a,b)
bandinds{DD,RR}(::ConcreteConversion{LaurentDirichlet{DD,RR},Laurent{DD,RR}}) = 0,2

function getindex{DD,RR}(C::ConcreteConversion{LaurentDirichlet{DD,RR},Laurent{DD,RR}},k::Integer,j::Integer)
    if k==1 && j==2
        one(eltype(C))
    elseif k==j || j==k+2
        one(eltype(C))
    else
        zero(eltype(C))
    end
end

conversion_rule{DD,RR}(b::LaurentDirichlet,a::Laurent{DD,RR})=b

differentiate{DD,RR}(f::Fun{LaurentDirichlet{DD,RR}}) = differentiate(Fun(f,Laurent))

for op in (:+,:-,:*)
    @eval begin
        $op{DD,RR}(f::Fun{Laurent{DD,RR}},g::Fun{LaurentDirichlet{DD,RR}}) = $op(f,Fun(g,Laurent))
        $op{DD,RR}(f::Fun{LaurentDirichlet{DD,RR}},g::Fun{Laurent{DD,RR}}) = $op(Fun(f,Laurent),g)
        $op{DD,RR}(f::Fun{LaurentDirichlet{DD,RR}},g::Fun{LaurentDirichlet{DD,RR}}) = $op(Fun(f,Laurent),Fun(g,Laurent))
    end
end


Base.real{DD,RR}(f::Fun{LaurentDirichlet{DD,RR}}) = real(Fun(f,Laurent))
Base.imag{DD,RR}(f::Fun{LaurentDirichlet{DD,RR}}) = imag(Fun(f,Laurent))


coefficients(v::AbstractVector,::Laurent,::LaurentDirichlet) = laurentdirichlettransform!(copy(v))
coefficients(v::AbstractVector,::LaurentDirichlet,::Laurent) = laurentidirichlettransform!(copy(v))

function laurentdirichlettransform!(w::AbstractVector)
    if length(w) > 1
        for k=length(w)-2:-1:1
            @inbounds w[k] -= w[k+2]
        end
        @inbounds w[1] -= w[2]
    end

    w
end

function laurentidirichlettransform!(w::AbstractVector)
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


struct CosDirichlet{DD,RR} <: Space{DD,RR}
    domain::DD
end

CosDirichlet(d::Domain) = CosDirichlet{typeof(d),real(eltype(d))}(d)
CosDirichlet() = CosDirichlet(PeriodicInterval())

spacescompatible(a::CosDirichlet,b::CosDirichlet) = domainscompatible(a,b)

canonicalspace(S::CosDirichlet)=CosSpace(domain(S))

Conversion(a::CosSpace,b::CosDirichlet) = ConcreteConversion(a,b)
bandinds{CS<:CosSpace,CD<:CosDirichlet}(::ConcreteConversion{CD,CS})=0,1
getindex{CS<:CosSpace,CD<:CosDirichlet}(C::ConcreteConversion{CD,CS},k::Integer,j::Integer) =
    (k==j||j==k+1)?one(eltype(C)):zero(eltype(C))


conversion_rule(b::CosDirichlet,a::CosSpace)=b
