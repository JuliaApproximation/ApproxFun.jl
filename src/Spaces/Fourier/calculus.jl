##Differentiation and integration

Base.sum{DD}(f::Fun{Laurent{DD}})=fouriersum(domain(f),f.coefficients)
Base.sum{DD}(f::Fun{Fourier{DD}})=fouriersum(domain(f),f.coefficients)

function linesum{DD}(f::Fun{Laurent{DD}})
    d=domain(f)
    if isa(d,Circle)
        sum(setcanonicaldomain(f))*d.radius
    else
        sum(f)
    end
end

linesum{DD<:Circle}(f::Fun{Fourier{DD}})=sum(setcanonicaldomain(f))*d.radius
linesum{DD<:PeriodicInterval}(f::Fun{Fourier{DD}})=sum(f) #TODO: Complex periodic interval


differentiate{DD<:PeriodicInterval}(f::Fun{Taylor{DD}}) = Fun(im*tocanonicalD(f,0)*taylor_diff(f.coefficients),f.space)
differentiate{DD<:PeriodicInterval}(f::Fun{Hardy{false,DD}}) = Fun(im*tocanonicalD(f,0)*hardyfalse_diff(f.coefficients),f.space)
differentiate{DD<:PeriodicInterval}(f::Fun{Laurent{DD}}) = Fun(im*tocanonicalD(f,0)*laurentdiff(f.coefficients),f.space)

differentiate{DD<:PeriodicInterval}(f::Fun{CosSpace{DD}}) = Fun(tocanonicalD(f,0)*cosspacediff(f.coefficients),SinSpace(domain(f)))
differentiate{DD<:PeriodicInterval}(f::Fun{SinSpace{DD}}) = Fun(tocanonicalD(f,0)*sinspacediff(f.coefficients),CosSpace(domain(f)))
differentiate{DD<:PeriodicInterval}(f::Fun{Fourier{DD}}) = Fun(tocanonicalD(f,0)*fourierdiff(f.coefficients),f.space)

differentiate{DD}(f::Fun{Laurent{DD}}) = Derivative(space(f))*f
differentiate{DD}(f::Fun{Fourier{DD}}) = Derivative(space(f))*f

function integrate{D}(f::Fun{Hardy{false,D}})
    if isa(domain(f),Circle) # drop -1 term if zero and try again
        @assert ncoefficients(f)==0 || abs(f.coefficients[1])<100eps()
        integrate(Fun(f,space(f)|(2:∞)))
    else  # Probably periodic itnerval
        Integral(space(f))*f
    end
end

function integrate{D}(f::Fun{Taylor{D}})
    if isa(domain(f),Circle)
        Integral(space(f))*f
    else  # Probably periodic itnerval  drop constant term if zero
        @assert ncoefficients(f)==0 || abs(f.coefficients[1])<100eps()
        Fun(integrate(Fun(f,space(f)|(2:∞))),space(f))
    end
end


function integrate{CS<:CosSpace}(f::Fun{CS})
    if isa(domain(f),Circle)
        error("Integrate not implemented for CosSpace on Circle")
    else  # Probably periodic itnerval, drop constant term if zero
        tol=1E-14 #TODO: smart tolerance.  Here relative is a bit tricky
                  # since this is called by Fourier integrate
        if abs(f.coefficients[1])<tol
            integrate(Fun(f,space(f)|(2:∞)))
        else
            d=domain(f)
            @assert isa(d,PeriodicInterval)
            x=Fun(identity,[first(d),last(d)])
            (f.coefficients[1]*x)⊕integrate(Fun(f,space(f)|(2:∞)))
        end
    end
end

function integrate{SS<:SinSpace}(f::Fun{SS})
    if isa(domain(f),Circle) # drop term containing z^(-1)
        integrate(Fun(f,space(f)|(2:∞)))
    else  # Probably periodic itnerval\
        Integral(space(f))*f
    end
end

#TODO: This is a hack to make sure Fourier maps to Fourier
# we don't have banded differentiate from CosSpace/SinSpace on a circle
for OP in (:differentiate,:integrate)
    @eval begin
        $OP{T,D<:PeriodicInterval}(f::Fun{Fourier{D},T})=$OP(f[2])⊕$OP(f[1])
        $OP{T,D<:Circle}(f::Fun{Fourier{D},T})=$OP(Fun(f,Laurent))
    end
end




fouriersum(d::PeriodicInterval,cfs)=cfs[1].*arclength(d)

function fouriersum{T}(d::Circle,cfs::Vector{T})
    if length(cfs)≥2
        2π*im*cfs[2]*d.radius
    else
        im*zero(T)
    end
end


# O(min(m,n)) Laurent line integral

function linebilinearform{T,D<:Circle}(f::Fun{Laurent{D},T},g::Fun{Laurent{D},T})
    @assert domain(f) == domain(g)
    u,v,mn = f.coefficients,g.coefficients,min(ncoefficients(f),ncoefficients(g))
    if mn > 1
        ret = u[1]*v[1]
        for i=2:2:mn-1
            ret += u[i]*v[i+1] + u[i+1]*v[i]
        end
        return arclength(domain(f))*ret
    elseif mn > 0
        return arclength(domain(f))*u[1]*v[1]
    else
        return zero(T)
    end
end

function bilinearform{T,D<:Circle}(f::Fun{Laurent{D},T},g::Fun{Laurent{D},T})
    @assert domain(f) == domain(g)
    u,v,mn = f.coefficients,g.coefficients,min(ncoefficients(f),ncoefficients(g))
    if mn > 2
        ret = u[1]*v[2] + u[2]*v[1]
        for i=3:2:mn-1
            ret += u[i]*v[i+1] + u[i+1]*v[i]
        end
        return complexlength(domain(f))*ret
    else
        return zero(T)
    end
end
