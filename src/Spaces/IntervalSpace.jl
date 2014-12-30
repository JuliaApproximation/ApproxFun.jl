
typealias IntervalSpace  RealSpace{Interval}     # We assume basis is real
canonicaldomain{T<:IntervalSpace}(::Type{T})=Interval()

Space{N<:Number}(d::Vector{N})=Space(Interval(d))

## Calculus

# the default domain space is higher to avoid negative ultraspherical spaces
Integral(d::IntervalDomain,n::Integer)=Integral(Ultraspherical{1}(d),n)

## Sigma

Σ(d::IntervalDomain)=Σ(JacobiWeight(-.5,-.5,Chebyshev(d)),Chebyshev(d))

function Σ(α::Number,β::Number,d::IntervalDomain)
    @assert α == β
    @assert int(α+.5) == α+.5
    @assert int(α+.5) >= 0
    Σ(JacobiWeight(α,β,Ultraspherical{int(α+.5)}(d)),Ultraspherical{int(α+.5)}(d))
end
Σ(α::Number,β::Number) = Σ(α,β,Interval())

## Evaluation

Evaluation(d::IntervalDomain,x::Union(Number,Bool),n...)=Evaluation(Chebyshev(d),x,n...)
Evaluation{T<:Number}(d::Vector{T},x::Union(Number,Bool),o::Integer)=Evaluation(Interval(d),x,o)
Evaluation(x::Union(Number,Bool))=Evaluation(Interval(),x,0)





function continuity{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T},order::Integer)

    m=length(d)
    B=zeros(Functional,m-1,m)
    
    for k=1:m-1
        B[k,k]=Evaluation(d[k],true,order)
        B[k,k+1]=-Evaluation(d[k+1],false,order)
    end
    B
end

function continuity{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T},kr::UnitRange)
    @assert first(kr)==0
    m=length(d)
    B=zeros(Functional,length(kr)*(m-1),m)
    for r in kr
        B[(m-1)*r+1:(m-1)*(r+1),:]=continuity(d,r)
    end
    B
end



function dirichlet{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2,m)
    B[1,1]=ldirichlet(d[1]);B[2,end]=rdirichlet(d[end])
    [B;
    continuity(d,0:1)]
end

function neumann{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2,m)
    B[1,1]=ldirichlet(d[1]);B[2,end]=rdirichlet(d[end])
    [B;
    continuity(d,0:1)]
end


function periodic{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2,m)
    B[1,1]=ldirichlet(d[1]);B[1,end]=-rdirichlet(d[end])
    B[2,1]=lneumann(d[1]);B[2,end]=-rneumann(d[end])    
    [B;
    continuity(d,0:1)]
end

function periodic{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2,m)
    B[1:2,1]=ivp(d[1])
    [B;
    continuity(d,0:1)]
end



## Orthogonal polynomials

abstract PolynomialSpace <: IntervalSpace

bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=domainspace(M)

