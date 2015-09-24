##Differentiation and integration

Base.sum{DD}(f::Fun{Laurent{DD}})=fouriersum(domain(f),f.coefficients)

function linesum{DD}(f::Fun{Laurent{DD}})
    d=domain(f)
    if isa(d,Circle)
        sum(Fun(f.coefficients,S(canonicaldomain(f))))*d.radius
    else
        sum(f)
    end
end

linesum{DD<:Circle}(f::Fun{Fourier{DD}})=sum(Fun(f.coefficients,Fourier(canonicaldomain(f))))*d.radius
linesum{DD<:PeriodicInterval}(f::Fun{Fourier{DD}})=sum(f) #TODO: Complex periodic interval



function integrate{D}(f::Fun{Hardy{false,D}})
    if isa(domain(f),Circle) # drop -1 term if zero and try again
        @assert length(f)==0 || abs(f.coefficients[1])<100eps()
        integrate(Fun(f,SliceSpace(space(f),1)))
    else  # Probably periodic itnerval
        Integral(space(f))*f
    end
end

function integrate{D}(f::Fun{Taylor{D}})
    if isa(domain(f),Circle)
        Integral(space(f))*f
    else  # Probably periodic itnerval  drop constant term if zero
        @assert length(f)==0 || abs(f.coefficients[1])<100eps()
        Fun(integrate(Fun(f,SliceSpace(space(f),1))),space(f))
    end
end


function integrate{CS<:CosSpace}(f::Fun{CS})
    if isa(domain(f),Circle)
        error("Integrate not implemented for CosSpace on Circle")
    else  # Probably periodic itnerval, drop constant term if zero
        tol=1E-14 #TODO: smart tolerance.  Here relative is a bit tricky
                  # since this is called by Fourier integrate
        if abs(f.coefficients[1])<tol
            integrate(Fun(f,SliceSpace(space(f),1)))
        else
            d=domain(f)
            @assert isa(d,PeriodicInterval)
            x=Fun(identity,[first(d),last(d)])
            (f.coefficients[1]*x)⊕integrate(Fun(f,SliceSpace(space(f),1)))
        end
    end
end

function integrate{SS<:SinSpace}(f::Fun{SS})
    if isa(domain(f),Circle) # drop term containing z^(-1)
        integrate(Fun(f,SliceSpace(space(f),1)))
    else  # Probably periodic itnerval\
        Integral(space(f))*f
    end
end

#TODO: This is a hack to make sure Fourier maps to Fourier
# we don't have banded differentiate from CosSpace/SinSpace on a circle
for OP in (:differentiate,:integrate)
    @eval begin
        $OP{T,D<:PeriodicInterval}(f::Fun{Fourier{D},T})=$OP(vec(f,2))⊕$OP(vec(f,1))
        $OP{T,D<:Circle}(f::Fun{Fourier{D},T})=$OP(Fun(f,Laurent))
    end
end




fouriersum(d::PeriodicInterval,cfs)=cfs[1].*length(d)

function fouriersum{T}(d::Circle,cfs::Vector{T})
    if length(cfs)≥2
        2π*im*cfs[2]*d.radius
    else
        im*zero(T)
    end
end
