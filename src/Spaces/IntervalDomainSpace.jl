
typealias IntervalDomainSpace  DomainSpace{Float64,Interval}     # We assume basis is real
canonicaldomain{T<:IntervalDomainSpace}(::Type{T})=Interval()

## Evaluation

function dirichlet{T<:Union(IntervalDomain,IntervalDomainSpace)}(d::Vector{T})
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

function neumann{T<:Union(IntervalDomain,IntervalDomainSpace)}(d::Vector{T})
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

function periodic{T<:Union(IntervalDomain,IntervalDomainSpace)}(d::Vector{T})
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


