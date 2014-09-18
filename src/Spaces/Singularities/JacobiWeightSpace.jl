

#Ultraspherical Spaces

export JacobiWeightSpace


immutable JacobiWeightSpace <: IntervalDomainSpace
    α::Float64
    β::Float64
    domain::Interval
end

JacobiWeightSpace(a::Number,b::Number,d)=JacobiWeightSpace(1.0a,1.0b,d)

spacescompatible(A::JacobiWeightSpace,B::JacobiWeightSpace)=A.α==B.α && A.β == B.β


jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
jacobiweight(sp::JacobiWeightSpace,x)=jacobiweight(sp.α,sp.β,tocanonical(sp,x))

evaluate(f::Fun{JacobiWeightSpace},x)=jacobiweight(space(f),x).*Fun(f.coefficients,domain(f))[x]

points(sp::JacobiWeightSpace,n)=fromcanonical(sp,chebyshevroots(n))
transform(sp::JacobiWeightSpace,vals::Vector)=chebyshevrootstransform(vals./jacobiweight(sp,points(sp,length(vals))))
itransform(sp::JacobiWeightSpace,cfs::Vector)=ichebyshevrootstransform(cfs).*jacobiweight(sp,points(sp,length(cfs)))



##TODO: paradigm for same space
spaceconversion(f::Vector,sp::JacobiWeightSpace,::ChebyshevSpace)=spaceconversion(f,sp,JacobiWeightSpace(0,0,domain(sp)))
spaceconversion(f::Vector,::ChebyshevSpace,sp::JacobiWeightSpace)=spaceconversion(f,JacobiWeightSpace(0,0,domain(sp)),sp)
function spaceconversion(f::Vector,sp1::JacobiWeightSpace,sp2::JacobiWeightSpace)
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β
    if c==α && d==β
        f
    elseif c>α && d>β
        spaceconversion(divide_singularity(f),JacobiWeightSpace(α+1,β+1,domain(sp1)),sp2)
    elseif c>α
        spaceconversion(divide_singularity(-1,f),JacobiWeightSpace(α+1,β,domain(sp1)),sp2)    
    elseif d>β
        spaceconversion(divide_singularity(1,f),JacobiWeightSpace(α,β+1,domain(sp1)),sp2)        
    else
        error("Need to implement decreasing jacobi")
    end
end

increase_jacobi_parameter(f)=Fun(f,JacobiWeightSpace(f.space.α+1,f.space.β+1,domain(f)))
increase_jacobi_parameter(s,f)=s==-1?Fun(f,JacobiWeightSpace(f.space.α+1,f.space.β,domain(f))):Fun(f,JacobiWeightSpace(f.space.α,f.space.β+1,domain(f)))


## Algebra

for op in (:/,:./)
    @eval begin
        ($op)(c::Number,f::Fun{JacobiWeightSpace})=Fun(($op)(c,Fun(f.coefficients)).coefficients,JacobiWeightSpace(-f.space.α,-f.space.β,domain(f)))        
    end
end

function .*{T,N}(f::Fun{JacobiWeightSpace,T},g::Fun{JacobiWeightSpace,N})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β    
    Fun((Fun(f.coefficients).*Fun(g.coefficients)).coefficients,JacobiWeightSpace(fα+gα,fβ+gβ,domain(f)))
end


function ./{T,N}(f::Fun{JacobiWeightSpace,T},g::Fun{JacobiWeightSpace,N})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β    
    Fun((Fun(f.coefficients)./Fun(g.coefficients)).coefficients,JacobiWeightSpace(fα-gα,fβ-gβ,domain(f)))
end

for op in (:.*,:./)
    @eval ($op){T,N,a}(f::Fun{UltrasphericalSpace{a},T},g::Fun{JacobiWeightSpace,N})=$op(Fun(f,JacobiWeightSpace(0,0,domain(f))),g)
    @eval ($op){T,N,a}(f::Fun{JacobiWeightSpace,N},g::Fun{UltrasphericalSpace{a},T})=$op(f,Fun(g,JacobiWeightSpace(0,0,domain(g)))) 
end


## Calculus

function Base.sum(f::Fun{JacobiWeightSpace})
    ##TODO: generalize
    α,β=f.space.α,f.space.β
    
    if α==β==.5
        fromcanonicalD(f,0.)*spaceconversion(f.coefficients,UltrasphericalSpace{1}(domain(f)))[1]*π/2
    elseif α==β==0.
        sum(Fun(f.coefficients,domain(f)))
    elseif α==β==-.5
        fromcanonicalD(f,0.)*π*f.coefficients[1]
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

