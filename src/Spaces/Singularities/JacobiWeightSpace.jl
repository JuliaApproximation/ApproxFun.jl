

#Ultraspherical Spaces

export JacobiWeightSpace


immutable JacobiWeightSpace{S<:IntervalDomainSpace} <: IntervalDomainSpace
    α::Float64
    β::Float64
    space::S
end

JacobiWeightSpace(a::Number,b::Number,d::IntervalDomainSpace)=JacobiWeightSpace(1.0a,1.0b,d)
JacobiWeightSpace(a::Number,b::Number,d::Domain)=JacobiWeightSpace(1.0a,1.0b,Space(d))
JacobiWeightSpace(a,b)=JacobiWeightSpace(a,b,ChebyshevSpace())

domain(S::JacobiWeightSpace)=domain(S.space)

spacescompatible(A::JacobiWeightSpace,B::JacobiWeightSpace)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)


jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
jacobiweight(sp::JacobiWeightSpace,x)=jacobiweight(sp.α,sp.β,tocanonical(sp,x))

evaluate{S}(f::Fun{JacobiWeightSpace{S}},x)=jacobiweight(space(f),x).*Fun(f.coefficients,space(f).space)[x]


## Use 1st kind points to avoid singularities
points(sp::JacobiWeightSpace{ChebyshevSpace},n)=fromcanonical(sp,chebyshevroots(n))
transform(sp::JacobiWeightSpace{ChebyshevSpace},vals::Vector)=chebyshevrootstransform(vals./jacobiweight(sp,points(sp,length(vals))))
itransform(sp::JacobiWeightSpace{ChebyshevSpace},cfs::Vector)=ichebyshevrootstransform(cfs).*jacobiweight(sp,points(sp,length(cfs)))

itransform(S::JacobiWeightSpace,cfs::Vector)=jacobiweight(S,points(S,length(cfs))).*itransform(S.space,cfs)


##TODO: paradigm for same space
spaceconversion(f::Vector,sp::JacobiWeightSpace,S2::ChebyshevSpace)=spaceconversion(f,sp,JacobiWeightSpace(0,0,S2))
spaceconversion(f::Vector,S2::ChebyshevSpace,sp::JacobiWeightSpace)=spaceconversion(f,JacobiWeightSpace(0,0,S2),sp)
function spaceconversion(f::Vector,sp1::JacobiWeightSpace{ChebyshevSpace},sp2::JacobiWeightSpace{ChebyshevSpace})
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β
    if c==α && d==β
        f
    elseif c>α && d>β
        spaceconversion(divide_singularity(f),JacobiWeightSpace(α+1,β+1,sp1.space),sp2)
    elseif c>α
        spaceconversion(divide_singularity(-1,f),JacobiWeightSpace(α+1,β,sp1.space),sp2)    
    elseif d>β
        spaceconversion(divide_singularity(1,f),JacobiWeightSpace(α,β+1,sp1.space),sp2)        
    else
        error("Need to implement decreasing jacobi")
    end
end

increase_jacobi_parameter(f)=Fun(f,JacobiWeightSpace(f.space.α+1,f.space.β+1,space(f).space))
increase_jacobi_parameter(s,f)=s==-1?Fun(f,JacobiWeightSpace(f.space.α+1,f.space.β,space(f).space)):Fun(f,JacobiWeightSpace(f.space.α,f.space.β+1,space(f).space))


## Algebra

for op in (:/,:./)
    @eval begin
        ($op){S}(c::Number,f::Fun{JacobiWeightSpace{S}})=Fun(($op)(c,Fun(f.coefficients)).coefficients,JacobiWeightSpace(-f.space.α,-f.space.β,space(f).space))        
    end
end

function .*{S,V}(f::Fun{JacobiWeightSpace{S}},g::Fun{JacobiWeightSpace{V}})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β    
    m=(Fun(f.coefficients,space(f).space).*Fun(g.coefficients,space(g).space))
    Fun(m.coefficients,JacobiWeightSpace(fα+gα,fβ+gβ,space(m)))
end


function ./{T,N}(f::Fun{JacobiWeightSpace{T}},g::Fun{JacobiWeightSpace{N}})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β    
    m=(Fun(f.coefficients,space(f).space)./Fun(g.coefficients,space(g).space))
    Fun(m.coefficients,JacobiWeightSpace(fα-gα,fβ-gβ,space(m)))
end

for op in (:.*,:./)
    ##TODO: Make general 
    @eval ($op){S}(f::Fun,g::Fun{JacobiWeightSpace{S}})=$op(Fun(f,JacobiWeightSpace(0,0,space(f))),g)
    @eval ($op){S}(f::Fun{JacobiWeightSpace{S}},g::Fun)=$op(f,Fun(g,JacobiWeightSpace(0,0,space(g)))) 
end


## Calculus

function Base.sum(f::Fun{JacobiWeightSpace{ChebyshevSpace}})
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




## Operators

function rangespace{J<:JacobiWeightSpace}(D::Derivative{J})
    S=D.space
    if S.α==S.β==0
       JacobiWeightSpace(0.,0.,rangespace(Derivative(S.space)))
    elseif S.α==0
       JacobiWeightSpace(0.,S.β-1,rangespace(Derivative(S.space)))
    elseif S.β==0
       JacobiWeightSpace(S.α-1,0.,rangespace(Derivative(S.space)))
    else
        #We assume the range is the same as the derivative
        # but really in general it should be
        # rangespace(S.α*(1-x) - S.β*(1-x) +(1-x.^2)*Derivative(S.space))
        # if multiplying by x changes space this needs to be redone
        JacobiWeightSpace(S.α-1,S.β-1,rangespace(Derivative(S.space)))
    end
end

function addentries!{J<:JacobiWeightSpace}(D::Derivative{J},A::ShiftArray,kr::Range)
    S=D.space
    if S.α==S.β==0
        addentries!(Derivative(S.space),A,kr)
    elseif S.α==0
        x=Fun(identity,S.space)
        DD=-S.β*I +(1-x)*Derivative(S.space)
        addentries!(DD,A,kr)    
    elseif S.β==0
        x=Fun(identity,S.space)
        DD=S.α*I +(1+x)*Derivative(S.space)
        addentries!(DD,A,kr)        
    else 
        x=Fun(identity,S.space)
        DD=S.α*(1-x) - S.β*(1+x) +(1-x.^2)*Derivative(S.space)
        addentries!(DD,A,kr)
    end
    
    A
end


