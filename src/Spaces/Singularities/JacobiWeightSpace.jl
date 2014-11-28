

#Ultraspherical Spaces

export JacobiWeightSpace


immutable JacobiWeightSpace{S<:IntervalDomainSpace} <: IntervalDomainSpace
    α::Float64
    β::Float64
    space::S
    function JacobiWeightSpace(α::Float64,β::Float64,space::S)
        if isa(space,JacobiWeightSpace)
            JacobiWeightSpace(α+space.α,β+space.β,space.space)
        else
            new(α,β,space)
        end
    end
end

JacobiWeightSpace{S<:IntervalDomainSpace}(a::Number,b::Number,d::S)=JacobiWeightSpace{S}(float64(a),float64(b),d)
JacobiWeightSpace(a::Number,b::Number,d::Domain)=JacobiWeightSpace(float64(a),float64(b),Space(d))
JacobiWeightSpace(a,b)=JacobiWeightSpace(a,b,ChebyshevSpace())
JacobiWeightSpace(a::Number,b::Number,d::Vector)=JacobiWeightSpace(float64(a),float64(b),Interval(d))


domain(S::JacobiWeightSpace)=domain(S.space)

spacescompatible(A::JacobiWeightSpace,B::JacobiWeightSpace)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)


jacobiweight(α,β,x)=(1+x).^α.*(1-x).^β
jacobiweight(sp::JacobiWeightSpace,x)=jacobiweight(sp.α,sp.β,tocanonical(sp,x))

function evaluate{S,T}(f::Fun{JacobiWeightSpace{S},T},x)
    tol=1.0E-14
    fv=Fun(f.coefficients,space(f).space)[x]
    if isa(fv,Number)&&abs(fv)<tol
        zero(T)
    else
        jacobiweight(space(f),x).*fv
    end
end


## Use 1st kind points to avoid singularities
points(sp::JacobiWeightSpace,n)=fromcanonical(sp,chebyshevroots(n))
transform(sp::JacobiWeightSpace,vals::Vector)=chebyshevrootstransform(vals./jacobiweight(sp,points(sp,length(vals))))
itransform(sp::JacobiWeightSpace,cfs::Vector)=ichebyshevrootstransform(cfs).*jacobiweight(sp,points(sp,length(cfs)))


plan_itransform(S::JacobiWeightSpace,n::Integer)=points(S,n)
itransform(S::JacobiWeightSpace,cfs::Vector)=itransform(S,cfs,plan_itransform(S,length(cfs)))
itransform(S::JacobiWeightSpace,cfs::Vector,pts::Vector)=jacobiweight(S,pts).*itransform(S.space,cfs,pts)

##TODO: paradigm for same space
function spaceconversion(f::Vector,sp1::JacobiWeightSpace,sp2::JacobiWeightSpace)
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β
    
    if isapprox(c,α) && isapprox(d,β)
        spaceconversion(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end
spaceconversion(f::Vector,sp::JacobiWeightSpace,S2::IntervalDomainSpace)=spaceconversion(f,sp,JacobiWeightSpace(0,0,S2))
spaceconversion(f::Vector,S2::IntervalDomainSpace,sp::JacobiWeightSpace)=spaceconversion(f,JacobiWeightSpace(0,0,S2),sp)
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



function canonicalspace(S::JacobiWeightSpace)
    if S.α==0 && S.β==0
        canonicalspace(S.space)
    else
        #TODO: promote singularities?
        JacobiWeightSpace(S.α,S.β,canonicalspace(S.space))
    end
end

## Algebra

for op in (:/,:./)
    @eval begin
        ($op){S}(c::Number,f::Fun{JacobiWeightSpace{S}})=Fun(($op)(c,Fun(f.coefficients)).coefficients,JacobiWeightSpace(-f.space.α,-f.space.β,space(f).space))        
    end
end

function .^{J<:JacobiWeightSpace}(f::Fun{J},k::Float64)
    S=space(f)
    g=Fun(coefficients(f),S.space)^k
    Fun(coefficients(g),JacobiWeightSpace(k*S.α,k*S.β,space(g)))
end

function .*{S,V}(f::Fun{JacobiWeightSpace{S}},g::Fun{JacobiWeightSpace{V}})
    @assert domainscompatible(f,g)
    fα,fβ=f.space.α,f.space.β
    gα,gβ=g.space.α,g.space.β    
    m=(Fun(f.coefficients,space(f).space).*Fun(g.coefficients,space(g).space))
    if isapprox(fα+gα,0)&&isapprox(fβ+gβ,0)
        m
    else
        Fun(m.coefficients,JacobiWeightSpace(fα+gα,fβ+gβ,space(m)))
    end
end


./{T,N}(f::Fun{JacobiWeightSpace{T}},g::Fun{JacobiWeightSpace{N}})=f*(1/g)

for op in (:.*,:./)
    ##TODO: Make general 
    @eval ($op){S}(f::Fun,g::Fun{JacobiWeightSpace{S}})=$op(Fun(f,JacobiWeightSpace(0,0,space(f))),g)
    @eval ($op){S}(f::Fun{JacobiWeightSpace{S}},g::Fun)=$op(f,Fun(g,JacobiWeightSpace(0,0,space(g)))) 
end

