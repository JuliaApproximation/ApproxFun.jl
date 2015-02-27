typealias LineSpace{T} IntervalMappedSpace{Chebyshev,Line{T}}
typealias PeriodicLineSpace{T} PeriodicMappedSpace{Fourier,PeriodicLine{T},RealBasis}
typealias PeriodicLineDirichlet{T} PeriodicMappedSpace{LaurentDirichlet,PeriodicLine{T},ComplexBasis}
typealias CurveSpace{S,T,DS} MappedSpace{S,Curve{S},T,DS}
typealias OpenCurveSpace{S} CurveSpace{S,RealBasis,Interval}
typealias ClosedCurveSpace{S,T} CurveSpace{S,T,PeriodicInterval}



# Transform form chebyshev U series to dirichlet-neumann U series
function uneumann_dirichlet_transform{T<:Number}(v::Vector{T})
    n=length(v)
    w=Array(T,n-4)

    for k = n-4:-1:1
        sc=(3+k)*(4+k)/((-2-k)*(1+k))
        w[k]=sc*v[k+4] 
        
        if k <= n-6
            w[k]-=sc*2*(4+k)/(5+k)*w[k+2]
        end
        if k <= n-8
            w[k]+=sc*((6+k)/(4+k))*w[k+4]
        end
    end
    
    w
end


# This takes a vector in dirichlet-neumann series on [-1,1]
# and return coefficients in T series that satisfy
# (1-x^2)^2 u' = f
function uneumannrange_xsqd{T<:Number}(v::Vector{T})
    n = length(v)
    w=Array(T,n+1)
    
    for k=n:-1:1
        sc=-((16*(1+k)*(2+k))/(k*(3+k)*(4+k)))
        w[k+1]=sc*v[k]
        
        if k <= n-2
            w[k+1]-=sc*(k*(4+k))/(8(k+1))*w[k+3]
        end
        
        if k <= n-4
            w[k+1]+=sc*((k*(k+4))/(16(k+2)))*w[k+5]
        end
    end
    w[1]=zero(T)
    
    w
end




#integration functions

integrate{LS,T}(f::Fun{LineSpace{LS},T})=linsolve([ldirichlet(),Derivative()],Any[0.,f];tolerance=length(f)^2*max(1,maximum(f.coefficients))*10E-13)
integrate{LS,RR<:Ray,TT,DS,T}(f::Fun{MappedSpace{LS,RR,TT,DS},T})=linsolve([BasisFunctional(1),Derivative(space(f))],Any[0.,f];tolerance=length(f)*10E-15)
function integrate{LS<:JacobiWeight,RR<:Ray,T,TT,DS}(f::Fun{MappedSpace{LS,RR,TT,DS},T})
    # x^k -> x^(k+1)  so +1,-1 to singularities 
    # if the last entry of f is close to zero wei use the same space
    D=Derivative(MappedSpace(space(f).domain,JacobiWeight(space(f).space.α+1,(space(f).space.β==0&&abs(last(f))≤1E-9)?0:(space(f).space.β-1),domain(space(f).space))))
    linsolve(D,f;tolerance=length(f)*1.0E-15)
end

function integrate{LS<:JacobiWeight,RR<:Line,T,TT,DS}(f::Fun{MappedSpace{LS,RR,TT,DS},T})
    # x^k -> x^(k+1)  so +1,-1 to singularities 
    # if the last entry of f is close to zero wei use the same space
    D=Derivative(MappedSpace(space(f).domain,JacobiWeight(space(f).space.α-1,space(f).space.β-1,domain(space(f).space))))
    linsolve([Evaluation(domainspace(D),0.),D],Any[0.,f];tolerance=length(f)*1.0E-15)
end


Base.cumsum{LS<:JacobiWeight,RR<:Ray,T,TT,DS}(f::Fun{MappedSpace{LS,RR,TT,DS},T})=integrate(f) # the choice of space is zero at 0


function Base.sum{LS,T}(f::Fun{LineSpace{LS},T})
    d=domain(f)
    if d.α==d.β==-.5
        p=Fun(f.coefficients,f.space.space)  # project to [-1,1]
        q=divide_singularity(p)              # divide by (1-x^2), result is in Chebyshev
        r=Fun(q.coefficients,JacobiWeight(-.5,-.5,Interval()))  # multiply by jacobi weight
        sum(r)
    else
        cf = integrate(f)
        last(cf) - first(cf)
    end
end



## Derivative

function Derivative{SS<:FunctionSpace,LD<:Line,T}(S::MappedSpace{SS,LD,T},order::Int)
    d=domain(S)
    @assert d.α==-1&&d.β==-1
    x=Fun(identity,S)    
    D1=Derivative(S.space)
    DS=SpaceOperator(D1,S,MappedSpace(domain(S),rangespace(D1)))

    M1=Multiplication(Fun(1,d),Space(d))
    Mx=Multiplication(x^2,Space(d))
    M1x=M1+Mx
    u=M1x\(2/π)  #tocanonicalD(S,x)=2/π*(1/(1+x^2))

    M=Multiplication(u,DS|>rangespace)

    D=DerivativeWrapper(M*DS,1)
    
    if order==1
        D
    else
        Derivative(rangespace(D),order-1)*D
    end    
end


function Derivative{SS<:FunctionSpace,RD<:Ray,T}(S::MappedSpace{SS,RD,T},order::Int)
    d=domain(S)
    @assert d.centre==0 && d.angle==0 && d.orientation
# x=Fun(identity,Ray())
# M=Multiplication(1+2x+x^2,Space(Ray()))
# u=M\2

    D1=Derivative(S.space)
    DS=SpaceOperator(D1,S,MappedSpace(domain(S),rangespace(D1)))
    
    u=Fun([0.75,-1.0,0.25],Ray())
    M=Multiplication(u,DS|>rangespace)

    D=DerivativeWrapper(M*DS,1)
    
    if order==1
        D
    else
        Derivative(rangespace(D),order-1)*D
    end    
end




function Derivative{SS<:FunctionSpace}(S::MappedSpace{SS,PeriodicLine{false}},order::Int)
    d=domain(S)
    @assert d.centre==0  && d.L==1.0
    
    a=Fun([1.,0,1],PeriodicInterval())
    M=Multiplication(a,space(a))
    DS=Derivative(space(a))
    D=SpaceOperator(M*DS,S,S)
    
    if order==1
        DerivativeWrapper(D,1)
    else
        DerivativeWrapper(Derivative(rangespace(D),order-1)*D,order)
    end    
end


