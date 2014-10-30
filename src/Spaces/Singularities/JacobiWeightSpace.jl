

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
    α,β=f.space.α,f.space.β    
    if α <= -1.0 && β <= -1.0
        fs = Fun(f.coefficients,f.space.space)
        d = domain(fs)
        return Inf*fromcanonicalD(f,0.)*(sign(fs[d.a])+sign(fs[d.b]))/2
    else
        n = length(f)
        c = zeros(n)
        c[1] = 2.^(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+2)
        if n > 1
            c[2] = c[1]*(α-β)/(α+β+2)
            for i=1:n-2
                c[i+2] = (2(α-β)*c[i+1]-(α+β-i+2)*c[i])/(α+β+i+2)
            end
        end
        return fromcanonicalD(f,0.)*dot(f.coefficients,c)
    end
end

function differentiate{J<:JacobiWeightSpace}(f::Fun{J})
    S=f.space
    d=domain(f)    
    ff=Fun(f.coefficients,S.space)    
    if S.α==S.β==0
        u=differentiate(ff)
        Fun(u.coefficients,JacobiWeightSpace(0.,0.,space(u)))
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)  
        u=-Mp*S.β*ff +(1-M).*differentiate(ff)
        Fun(u.coefficients,JacobiWeightSpace(0.,S.β-1,space(u)))
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)  
        u=Mp*S.α*ff +(1+M).*differentiate(ff)
        Fun(u.coefficients,JacobiWeightSpace(S.α-1,0.,space(u)))        
    else 
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a) 
        u=(Mp*S.α)*(1-M).*ff- (Mp*S.β)*(1+M).*ff +(1-M.^2).*differentiate(ff)
        Fun(u.coefficients,JacobiWeightSpace(S.α-1,S.β-1,space(u)))        
    end
end




## Operators


function Derivative(S::JacobiWeightSpace)
    d=domain(S)
    @assert isa(d,Interval)

    if S.α==S.β==0
        DerivativeWrapper(SpaceOperator(Derivative(S.space),S,JacobiWeightSpace(0.,0.,rangespace(Derivative(S.space)))),1)
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)            
        DD=(-Mp*S.β)*I +(1-M)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(0.,S.β-1,rangespace(DD))),1)
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)        
        DD=(Mp*S.α)*I +(1+M)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(S.α-1,0.,rangespace(DD))),1)
    else 
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        DD=(Mp*S.α)*(1-M) - (Mp*S.β)*(1+M) +(1-M.^2)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(S.α-1,S.β-1,rangespace(DD))),1)
    end

end

function Derivative(S::JacobiWeightSpace,k::Integer)
    if k==1
        Derivative(S)
    else
        D=Derivative(S)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end




## Multiplication

addentries!{T,S<:JacobiWeightSpace}(M::Multiplication{T,S},A::ShiftArray,kr::Range)=addentries!(Multiplication(M.f,domainspace(M).space),A,kr)

addentries!{T<:JacobiWeightSpace,S<:JacobiWeightSpace}(M::Multiplication{T,S},A::ShiftArray,kr::Range)=addentries!(Multiplication(Fun(M.f.coefficients,space(M.f).space),domainspace(M).space),A,kr)
rangespace{T<:JacobiWeightSpace,S<:JacobiWeightSpace}(M::Multiplication{T,S})=JacobiWeightSpace(space(M.f).α+M.space.α,space(M.f).β+M.space.β,rangespace(Multiplication(M.f,domainspace(M).space)))



## Conversion

maxspace(A::JacobiWeightSpace,B::JacobiWeightSpace)=JacobiWeightSpace(min(A.α,B.α),min(A.β,B.β),maxspace(A.space,B.space))
minspace(A::JacobiWeightSpace,B::JacobiWeightSpace)=JacobiWeightSpace(max(A.α,B.α),max(A.β,B.β),minspace(A.space,B.space))
maxspace(A::IntervalDomainSpace,B::JacobiWeightSpace)=maxspace(JacobiWeightSpace(0.,0.,A),B)
maxspace(A::JacobiWeightSpace,B::IntervalDomainSpace)=maxspace(A,JacobiWeightSpace(0.,0.,B))

isapproxinteger(x)=isapprox(x,int(x))



function Conversion(A::JacobiWeightSpace,B::JacobiWeightSpace)
    @assert isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)
    
    if A.space==B.space
        d=domain(A)        
        x=Fun(identity,d)
        M=tocanonical(d,x)
        m=(1+M).^int(A.α-B.α).*(1-M).^int(A.β-B.β)
        SpaceOperator(Multiplication(m,B.space),A,B)# Wrap the operator with the correct spaces
    elseif isapprox(A.α,B.α) && isapprox(A.β,B.β)
        SpaceOperator(Conversion(A.space,B.space),A,B)
    else
        d=domain(A)        
        x=Fun(identity,d)
        M=tocanonical(d,x)    
        C=Conversion(A.space,B.space)
        m=(1+M).^int(A.α-B.α).*(1-M).^int(A.β-B.β)
        SpaceOperator(Multiplication(m,B.space)*C,A,B)
    end        
end


isapproxleq(a,b)=(a<=b || isapprox(a,b))
# return the space that has banded Conversion to the other, or NoSpace
function conversion_rule(A::JacobiWeightSpace,B::JacobiWeightSpace)
    if isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)    
        ct=conversion_type(A.space,B.space)
        if ct == B.space && isapproxleq(A.α,B.α) && isapproxleq(A.β,B.β)
            return B
        elseif ct == A.space && isapproxleq(B.α,A.α) && isapproxleq(B.β,A.β)
            return A
        end
    end

    return NoSpace()
end


## Evaluation

function  Base.getindex{J<:JacobiWeightSpace}(op::Evaluation{J,Bool},kr::Range)
    S=op.space
    @assert op.order<=1
    d=domain(op)
    @assert isa(d,Interval)
    if op.x
        @assert S.β>=0
        if S.β==0
            if op.order==0
                2^S.α*getindex(Evaluation(S.space,op.x),kr)
            else #op.order ===1
                2^S.α*getindex(Evaluation(S.space,op.x,1),kr)+(tocanonicalD(d,d.a)*S.α*2^(S.α-1))*getindex(Evaluation(S.space,op.x),kr)
            end
        else
            @assert op.order==0
            zeros(kr)
        end
    else
        @assert S.α>=0
        if S.α==0
            if op.order==0
                2^S.β*getindex(Evaluation(S.space,op.x),kr)
            else #op.order ===1
                2^S.β*getindex(Evaluation(S.space,op.x,1),kr)-(tocanonicalD(d,d.a)*S.β*2^(S.β-1))*getindex(Evaluation(S.space,op.x),kr)
            end
        else
            @assert op.order==0        
            zeros(kr)
        end    
    end
end

