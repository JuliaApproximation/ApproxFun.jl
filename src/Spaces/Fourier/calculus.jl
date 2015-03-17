##Differentiation and integration

export linesum


Base.sum(f::Fun{Laurent})=fouriersum(domain(f),f.coefficients)
function linesum{S<:PeriodicSpace}(f::Fun{S})
    d=domain(f)
    if isa(d,Circle)
        sum(Fun(f.coefficients,S(canonicaldomain(f))))*d.radius
    else
        sum(f)
    end
end

function integrate(f::Fun{Hardy{false}})
    if isa(domain(f),Circle) # drop -1 term if zero and try again
        @assert length(f)==0 || abs(f.coefficients[1])<100eps()
        integrate(Fun(f,SliceSpace(space(f),1)))
    else  # Probably periodic itnerval
        Integral(space(f))*f
    end
end

function integrate(f::Fun{Taylor})
    if isa(domain(f),Circle)
        Integral(space(f))*f
    else  # Probably periodic itnerval  drop constant term if zero
        @assert length(f)==0 || abs(f.coefficients[1])<100eps()
        Fun(integrate(Fun(f,SliceSpace(space(f),1))),space(f))
    end
end


function integrate(f::Fun{CosSpace})
    if isa(domain(f),Circle)
        error("Integrate not implemented for CosSpace on Circle")
    else  # Probably periodic itnerval, drop constant term if zero
        if isapprox(f.coefficients[1],0)
            integrate(Fun(f,SliceSpace(space(f),1)))
        else
            d=domain(f)
            @assert isa(d,PeriodicInterval)
            x=Fun(identity,[first(d),last(d)])
            (f.coefficients[1]*x)⊕integrate(Fun(f,SliceSpace(space(f),1)))
        end
    end
end

function integrate(f::Fun{SinSpace})
    if isa(domain(f),Circle) # drop term containing z^(-1)
        integrate(Fun(f,SliceSpace(space(f),1)))
    else  # Probably periodic itnerval\
        Integral(space(f))*f
    end
end

#TODO: This is a hack to make sure Fourier maps to Fourier
# we don't have banded differentiate from CosSpace/SinSpace on a circle
for OP in (:differentiate,:integrate)
    @eval $OP{T}(f::Fun{Fourier,T})=isa(domain(f),PeriodicInterval)?($OP(vec(f,2))⊕$OP(vec(f,1))):$OP(Fun(f,Laurent))
end




fouriersum(d::PeriodicInterval,cfs)=cfs[1].*length(d)

function fouriersum{T}(d::Circle,cfs::Vector{T})
    if length(cfs)≥2
        2π*im*cfs[2]*d.radius
    else
        im*zero(T)
    end
end

