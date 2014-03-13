function givensreduceab!{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer)
    a=datagetindex(B,k1,j1)
    b=datagetindex(B,k2,j1)
    
    if b == 0.
        return one(T),zero(T)
    end    
    

    
    sq=sqrt(abs2(a) + abs2(b))    
    a=a/sq;b=b/sq
    
    #TODO: Assuming that left rows are already zero
    
    ir1=indexrange(B,k1)::Range1{Int64}
    ir2=indexrange(B,k2)::Range1{Int64}    
    
    for j = j1:ir1[end]
        B1 = datagetindex(B,k1,j)
        B2 = datagetindex(B,k2,j)
        
        B[k1,j],B[k2,j]= a*B1 + b*B2,-b*B1 + a*B2
    end
    
    for j=ir1[end]+1:ir2[end]
        B1 = fillgetindex(B,k1,j)
        B2 = datagetindex(B,k2,j)
        
        B[k2,j]=a*B2 - b*B1
    end
    
    for j=1:numbcs(B)
        B1 = getfilldata(B,k1,j)
        B2 = getfilldata(B,k2,j)
    
        setfilldata!(B, a*B1 + b*B2,k1,j)
        setfilldata!(B,-b*B1 + a*B2,k2,j)    
    end
    

    a::T,b::T
end


#assume dimension of kernel is range space
Base.null(A::BandedOperator)=null(A,rangespace(A))


function applygivens!(Q,k,a,b)
    for j=1:length(Q[1])
        Q[1][j],Q[k][j]=a*Q[1][j]+b*Q[k][j],-b*Q[1][j]+a*Q[k][j]
    end
end


#d is number of elements in the kernel
function Base.null{T<:Number}(A::BandedOperator{T},d,maxit=Inf)
    M=MutableAlmostBandedOperator([A'])
    n=100  
    resizedata!(M,n)
    m=bandrange(A)[end]
    Q=Vector{T}[zeros(T,n) for j=1:m]
    for j=1:m
        Q[j][j]=one(T)
    end
    
    ret=Q[1]
     
    k=0
    
    while norm(M.data.data[k+1:k+d,:])>eps()
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
    
    
    IFun[Fun(Q[j][1:k],domain(A)) for j=1:d]
end

