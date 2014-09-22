

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
spaceconversion(f::Vector,sp::JacobiWeightSpace,S2::JacobiWeightSpace)=error("Still need to implement")
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

function Derivative{JS<:JacobiWeightSpace}(J::JS,k::Integer)
    D=Derivative{JS,Float64}(J,1)
    if k==1
        D
    else
        TimesOperator(Derivative(rangespace(D),k-1),D)
    end
end


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
    @assert D.order ==1
    
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


## Multiplication

addentries!{T,S<:JacobiWeightSpace}(M::Multiplication{T,S},A::ShiftArray,kr::Range)=addentries!(Multiplication(M.f,domainspace(M).space),A,kr)

addentries!{T<:JacobiWeightSpace,S<:JacobiWeightSpace}(M::Multiplication{T,S},A::ShiftArray,kr::Range)=addentries!(Multiplication(Fun(M.f.coefficients,space(M.f).space),domainspace(M).space),A,kr)
rangespace{T<:JacobiWeightSpace,S<:JacobiWeightSpace}(M::Multiplication{T,S})=JacobiWeightSpace(space(M.f).α+M.space.α,space(M.f).β+M.space.β,rangespace(Multiplication(M.f,domainspace(M).space)))



## Conversion

maxspace(A::JacobiWeightSpace,B::JacobiWeightSpace)=JacobiWeightSpace(min(A.α,B.α),min(A.β,B.β),maxspace(A.space,B.space))
minspace(A::JacobiWeightSpace,B::JacobiWeightSpace)=JacobiWeightSpace(max(A.α,B.α),max(A.β,B.β),minspace(A.space,B.space))

function addentries!{Y<:JacobiWeightSpace,W<:JacobiWeightSpace}(C::Conversion{Y,W},SA::ShiftArray,kr::Range)
    @assert domain(C)==Interval()  ##TODO: General domains
    A=C.domainspace;B=C.rangespace
    @assert isinteger(A.α-B.α) && isinteger(A.β-B.β)
    if A.space==B.space
        x=Fun(identity)
        m=(1+x).^int(A.α-B.α).*(1-x).^int(A.β-B.β)
        addentries!(Multiplication(m,B.space),SA,kr)
    elseif A.α==B.α && A.β==B.β
        addentries!(Conversion(A.space,B.space),SA,kr)
    else
        C=Conversion(A.space,B.space)
        x=Fun(identity)
        m=(1+x).^int(A.α-B.α).*(1-x).^int(A.β-B.β)
        addentries!(Multiplication(m,B.space)*C,SA,kr)            
    end
end

function bandinds{Y<:JacobiWeightSpace,W<:JacobiWeightSpace}(C::Conversion{Y,W})
    A=C.domainspace;B=C.rangespace
    @assert isinteger(A.α-B.α) && isinteger(A.β-B.β)
    if A.space==B.space
        x=Fun(identity)
        m=(1+x).^int(A.α-B.α).*(1-x).^int(A.β-B.β)
        bandinds(Multiplication(m,B.space))
    elseif A.α==B.α && A.β==B.β
        bandinds(Conversion(A.space,B.space))
    else
        C=Conversion(A.space,B.space)
        x=Fun(identity)
        m=(1+x).^int(A.α-B.α).*(1-x).^int(A.β-B.β)
        bandinds(Multiplication(m,B.space)*C)
    end
end


## Evaluation

function  Base.getindex{J<:JacobiWeightSpace}(op::Evaluation{J,Bool},kr::Range)
    S=op.space
    if op.x
        @assert S.β>=0
        if S.β==0
            getindex(Evaluation(S.space,op.x),kr)
        else
            zeros(kr)
        end
    else
        @assert S.α>=0
        if S.α==0
            getindex(Evaluation(S.space,op.x),kr)
        else
            zeros(kr)
        end    
    end
end

