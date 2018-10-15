##Differentiation and integration


Base.sum(f::Fun{Laurent{DD,RR}}) where {DD<:PeriodicSegment,RR} = coefficient(f,1).*arclength(domain(f))
Base.sum(f::Fun{Laurent{DD,RR}}) where {DD<:Circle,RR} = coefficient(f,3 - domain(f).orientation).*complexlength(domain(f))


Base.sum(f::Fun{Fourier{DD,RR}}) where {DD<:PeriodicSegment,RR} = coefficient(f,1).*arclength(domain(f))
Base.sum(f::Fun{Fourier{DD,RR}}) where {DD<:Circle,RR} =
    (im*coefficient(f,2) + coefficient(f,3))/2*complexlength(domain(f))


function linesum(f::Fun{Laurent{DD,RR}}) where {DD,RR}
    d=domain(f)
    if isa(d,Circle)
        sum(setcanonicaldomain(f))*d.radius
    else
        sum(f)
    end
end

linesum(f::Fun{Fourier{DD,RR}}) where {DD<:Circle,RR} = sum(setcanonicaldomain(f))*domain(f).radius
linesum(f::Fun{Fourier{DD,RR}}) where {DD<:PeriodicSegment,RR} = sum(f) #TODO: Complex periodic interval


differentiate(f::Fun{Taylor{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(f.space,im*tocanonicalD(f,0)*taylor_diff(f.coefficients))
differentiate(f::Fun{Hardy{false,DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(f.space,im*tocanonicalD(f,0)*hardyfalse_diff(f.coefficients))
differentiate(f::Fun{Laurent{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(f.space,im*tocanonicalD(f,0)*laurentdiff(f.coefficients))

differentiate(f::Fun{CosSpace{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(SinSpace(domain(f)),tocanonicalD(f,0)*cosspacediff(f.coefficients))
differentiate(f::Fun{SinSpace{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(CosSpace(domain(f)),tocanonicalD(f,0)*sinspacediff(f.coefficients))
differentiate(f::Fun{Fourier{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    Fun(f.space,tocanonicalD(f,0)*fourierdiff(f.coefficients))

differentiate(f::Fun{Laurent{DD,RR}}) where {DD,RR} = Derivative(space(f))*f
differentiate(f::Fun{Fourier{DD,RR}}) where {DD,RR} = Derivative(space(f))*f

function integrate(f::Fun{Hardy{false,D,R}}) where {D,R}
    if isa(domain(f),Circle) # drop -1 term if zero and try again
        @assert ncoefficients(f)==0 || abs(f.coefficients[1])<100eps()
        integrate(Fun(f,space(f)|(2:∞)))
    else  # Probably periodic itnerval
        Integral(space(f))*f
    end
end

function integrate(f::Fun{Taylor{D,R}}) where {D,R}
    if isa(domain(f),Circle)
        Integral(space(f))*f
    else  # Probably periodic itnerval  drop constant term if zero
        @assert ncoefficients(f)==0 || abs(f.coefficients[1])<100eps()
        Fun(integrate(Fun(f,space(f)|(2:∞))),space(f))
    end
end


Base.sum(f::Fun{CosSpace{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    f.coefficients[1]*complexlength(domain(f))

linesum(f::Fun{CosSpace{DD,RR}}) where {DD<:PeriodicSegment,RR} =
    f.coefficients[1]*arclength(domain(f))



function integrate(f::Fun{CS}) where CS<:CosSpace
    if isa(domain(f),Circle)
        error("Integrate not implemented for CosSpace on Circle")
    else  # Probably periodic itnerval, drop constant term if zero
        tol=1E-14 #TODO: smart tolerance.  Here relative is a bit tricky
                  # since this is called by Fourier integrate
        if abs(f.coefficients[1])<tol
            integrate(Fun(f,space(f)|(2:∞)))
        else
            d=domain(f)
            @assert isa(d,PeriodicSegment)
            x=Fun(identity, Interval(d))
            (f.coefficients[1]*x)⊕integrate(Fun(f,space(f)|(2:∞)))
        end
    end
end

function integrate(f::Fun{SS}) where SS<:SinSpace
    if isa(domain(f),Circle) # drop term containing z^(-1)
        integrate(Fun(f,space(f)|(2:∞)))
    else  # Probably periodic itnerval\
        Integral(space(f))*f
    end
end

#TODO: This is a hack to make sure Fourier maps to Fourier
# we don't have banded differentiate from CosSpace/SinSpace on a circle
for OP in (:differentiate,:integrate)
    @eval $OP(f::Fun{Fourier{D,R},T}) where {T,D<:Circle,R} = $OP(Fun(f,Laurent))
end

integrate(f::Fun{Fourier{D,R},T}) where {T,D<:PeriodicSegment,R} =
    integrate(component(f,2))⊕integrate(component(f,1))




# O(min(m,n)) Laurent line integral

function linebilinearform(f::Fun{Laurent{D,R},T},g::Fun{Laurent{D,R},T}) where {T,D<:Circle,R}
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

function bilinearform(f::Fun{Laurent{D,R},T},g::Fun{Laurent{D,R},T}) where {T,D<:Circle,R}
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
