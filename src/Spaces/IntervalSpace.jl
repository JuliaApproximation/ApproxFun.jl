
typealias IntervalSpace  RealSpace{Interval}     # We assume basis is real
canonicaldomain{T<:IntervalSpace}(::Type{T})=Interval()

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


function dirichlet{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2m,m)
    B[1,1]=dirichlet(d[1])[1]
    B[2,end]=dirichlet(d[end])[end]
    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end

function neumann{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2m,m)
    B[1,1]=neumann(d[1])[1]
    B[2,end]=neumann(d[end])[end]
    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end

function periodic{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2m,m)
    B[1,1]=dirichlet(d[1])[1]
    B[1,end]=-dirichlet(d[end])[end]
    B[2,1]=neumann(d[1])[1]
    B[2,end]=-neumann(d[end])[end]

    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end



function ivp{T<:Union(IntervalDomain,IntervalSpace)}(d::Vector{T})
    m=length(d)
    B=zeros(Functional,2m,m)
    B[1:2,1]=ivp(d[1])
    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end
