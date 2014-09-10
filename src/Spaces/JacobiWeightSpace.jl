

#Ultraspherical Spaces



immutable JacobiWeightSpace{α,β} <: IntervalDomainSpace
    domain::Interval
end

jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
jacobiweight{α,β}(sp::JacobiSpace{α,β},x)=jacobiweight(α,β,tocanonical(sp,x))

evaluate{T,J<:JacobiWeightSpace}(f::IFun{T,J},x)=jacobiweight(space(f),x).*IFun(f.coefficients,domain(f))[x]

#TODO: transform and use first kind points
itransform{J<:JacobiWeightSpace}(sp::J,cfs::Vector)=itransform(ChebyshevSpace(domain(sp)),cfs).*jacobiweight(sp,points(sp,length(cfs)))



##TODO: paradigm for same space
spaceconversion(f::Vector,::JacobiWeightSpace{0,0},::ChebyshevSpace)=f
spaceconversion(f::Vector,::ChebyshevSpace,::JacobiWeightSpace{0,0})=f

function Base.sum{T,α,β}(f::IFun{T,JacobiWeightSpace{α,β}})
    ##TODO: generalize

    if α==β==.5
        fromcanonicalD(f,0.)*coefficients(f.fun,UltrasphericalSpace{1}(domain(f)))[1]*π/2
    elseif α==β==0.
        sum(Fun(f.coefficients,domain(f)))
    elseif α==β==-.5
        fromcanonicalD(f,0.)*π*fun.coefficients[1]
    elseif α<0. && β<0.
        #TODO: should be < -1.
        sum(increase_jacobi_parameter(f))
    elseif α < 0
        sum(increase_jacobi_parameter(-1,f))
    elseif  β < 0
        sum(increase_jacobi_parameter(+1,f))    
    else
        error("sum not implemented for all Jacobi parameters")
    end
end