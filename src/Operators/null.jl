


#assume dimension of kernel is range space
Base.null(A::BandedOperator)=null(A,order(rangespace(A)))


function applygivens!(Q,k,a,b)
    for j=1:length(Q[1])
        Q[1][j],Q[k][j]=a*Q[1][j]+b*Q[k][j],-b*Q[1][j]+a*Q[k][j]
    end
end


#d is number of elements in the kernel
function Base.null{T<:Number}(A::BandedOperator{T},d,maxit=Inf)
    M=MutableAlmostBandedOperator([A'])
    m=bandinds(A)[end]
    n=m+100  
    resizedata!(M,n)
    Q=Vector{T}[zeros(T,n) for j=1:m]
    for j=1:m
        Q[j][j]=one(T)
    end
    
    ret=Q[1]
     
    k=0
    
    while norm(M.data.data[k+1:k+d,:])>eps()  && k <= maxit
        k+=1
        
        if k+m >= n
            n=2n
            resizedata!(M,2n)
            
            for j=1:m
                Q[j]=[Q[j],zeros(T,n)]
            end
        end
        
        for j=1:m-1
            applygivens!(Q,1+j,givensreduceab!(M,k,k+j,k)... )
        end
        
        ret=Q[1]
            
        for j=1:m-1
            Q[j]=Q[j+1]
        end
        
        a,b=givensreduceab!(M,k,k+m,k)   
        
        Q[m] = -b*ret
        Q[m][k+m] += a
    end
    
    
    Fun[Fun(Q[j][1:k],domain(A)) for j=1:d]
end


function Base.null{T<:BandedOperator}(A::Array{T,2},d)
    ret = null(interlace(A),d)
    
    [Fun(ret[1].coefficients[1:2:end],ret[1].domain) Fun(ret[2].coefficients[1:2:end],ret[1].domain);
    Fun(ret[1].coefficients[2:2:end],ret[1].domain) Fun(ret[2].coefficients[2:2:end],ret[1].domain)]
end

function Base.null{T<:BandedOperator}(A::Array{T,2})
    d = 0
    for k=1:size(A,1)
        d += max(rangespace(A[k,1]),rangespace(A[k,2]))
    end
    
    null(A,d)
end
