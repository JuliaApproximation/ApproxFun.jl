export Evaluation

## Evaluation constructors

# M = Bool if endpoint
type Evaluation{S<:FunctionSpace,M<:Union(Number,Bool),T<:Number} <: Functional{T}
    space::S
    x::M
    order::Int
end
Evaluation(sp::AnySpace,x::Bool)=Evaluation{AnySpace,Bool,Float64}(sp,x,0)
Evaluation(sp::AnySpace,x::Bool,k::Integer)=Evaluation{AnySpace,Bool,Float64}(sp,x,k)
Evaluation{M,S<:IntervalDomainSpace}(sp::S,x::M,order::Integer)=Evaluation{S,M,Float64}(sp,x,order)
Evaluation{M,S<:PeriodicDomainSpace}(sp::S,x::M,order::Integer)=Evaluation{S,M,Complex{Float64}}(sp,x,order)


Evaluation(d::FunctionSpace,x::Union(Number,Bool))=Evaluation(d,x,0)
Evaluation(d::IntervalDomain,x::Union(Number,Bool),n...)=Evaluation(ChebyshevSpace(d),x,n...)
Evaluation(d::PeriodicDomain,x::Number,n...)=Evaluation(LaurentSpace(d),complex(x),n...)
Evaluation{T<:Number}(d::Vector{T},x::Union(Number,Bool),o::Integer)=Evaluation(Interval(d),x,o)
Evaluation(x::Union(Number,Bool))=Evaluation(Interval(),x,0)

domainspace(E::Evaluation)=E.space
domain(E::Evaluation)=domain(E.space)

promotedomainspace(E::Evaluation,sp::FunctionSpace)=Evaluation(sp,E.x,E.order)


## Convenience routines


evaluate(d::IntervalDomain,x)=Evaluation(d,x)
ldirichlet(d::IntervalDomain)=Evaluation(d,false)
rdirichlet(d::IntervalDomain)=Evaluation(d,true)
lneumann(d::IntervalDomain)=Evaluation(d,false,1)
rneumann(d::IntervalDomain)=Evaluation(d,true,1)


ldirichlet(d::FunctionSpace)=Evaluation(d,false)
rdirichlet(d::FunctionSpace)=Evaluation(d,true)
lneumann(d::FunctionSpace)=Evaluation(d,false,1)
rneumann(d::FunctionSpace)=Evaluation(d,true,1)


ldiffbc(d::IntervalDomain,k::Integer) = Evaluation(d,false,k)
rdiffbc(d::IntervalDomain,k::Integer) = Evaluation(d,true,k)
ldiffbc(d::FunctionSpace,k::Integer) = Evaluation(d,false,k)
rdiffbc(d::FunctionSpace,k::Integer) = Evaluation(d,true,k)


dirichlet(d::Union(IntervalDomain,FunctionSpace))=[ldirichlet(d),rdirichlet(d)]
neumann(d::Union(IntervalDomain,FunctionSpace))=[lneumann(d),rneumann(d)]
diffbcs(d::Union(IntervalDomain,FunctionSpace),k::Integer) = [ldiffbc(d,k),rdiffbc(d,k)]


for op in (:rdirichlet,:ldirichlet,:dirichlet,:lneumann,:rneumann,:neumann)
    @eval $op()=$op(AnySpace())
end

for op in (:ldiffbc,:rdiffbc,:diffbcs)
    @eval $op(k::Integer)=$op(AnySpace(),k)
end

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