export EvaluationFunctional

## EvaluationFunctional constructors

# M = Bool if endpoint
type EvaluationFunctional{T<:Number,M<:Union(Number,Bool),S<:FunctionSpace} <: Functional{T}
    space::S
    x::M
    order::Int
end

EvaluationFunctional{M,S<:IntervalDomainSpace}(sp::S,x::M,order::Integer)=EvaluationFunctional{Float64,M,S}(sp,x,order)
EvaluationFunctional{M,S<:PeriodicDomainSpace}(sp::S,x::M,order::Integer)=EvaluationFunctional{Complex{Float64},M,S}(sp,x,order)


EvaluationFunctional(d::FunctionSpace,x::Union(Number,Bool))=EvaluationFunctional(d,x,0)
EvaluationFunctional(d::IntervalDomain,x::Union(Number,Bool),n...)=EvaluationFunctional(ChebyshevSpace(d),x,n...)
EvaluationFunctional(d::PeriodicDomain,x::Number,n...)=EvaluationFunctional(LaurentSpace(d),complex(x),n...)
EvaluationFunctional{T<:Number}(d::Vector{T},x::Union(Number,Bool),o::Integer)=EvaluationFunctional(Interval(d),x,o)
EvaluationFunctional(x::Union(Number,Bool))=EvaluationFunctional(Interval(),x,0)

domainspace(E::EvaluationFunctional)=E.space
domain(E::EvaluationFunctional)=domain(E.space)


## Convenience routines


evaluate(d::IntervalDomain,x)=EvaluationFunctional(d,x)
ldirichlet(d::IntervalDomain)=EvaluationFunctional(d,false)
rdirichlet(d::IntervalDomain)=EvaluationFunctional(d,true)
lneumann(d::IntervalDomain)=EvaluationFunctional(d,false,1)
rneumann(d::IntervalDomain)=EvaluationFunctional(d,true,1)


ldirichlet(d::IntervalDomainSpace)=EvaluationFunctional(d,false)
rdirichlet(d::IntervalDomainSpace)=EvaluationFunctional(d,true)
lneumann(d::IntervalDomainSpace)=EvaluationFunctional(d,false,1)
rneumann(d::IntervalDomainSpace)=EvaluationFunctional(d,true,1)


dirichlet(d::Union(IntervalDomain,IntervalDomainSpace))=[ldirichlet(d),rdirichlet(d)]
neumann(d::Union(IntervalDomain,IntervalDomainSpace))=[lneumann(d),rneumann(d)]


function dirichlet{T<:IntervalDomain}(d::Vector{T})
    m=length(d)
    B=zeros(Operator,2m,m)
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

function periodic{T<:IntervalDomain}(d::Vector{T})
    m=length(d)
    B=zeros(Operator,2m,m)
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