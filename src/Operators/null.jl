function givensreduceab!{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer,n)
    a=datagetindex(B,k1,j1)
    b=datagetindex(B,k2,j1)
    
    if b == 0.
        return B,1.,0.
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
    

    B,a,b
end




function Base.null(A::BandedOperator)
    M=MutableAlmostBandedOperator([A'])
    n=200  #TODO: change 200
    resizedata!(M,n)
    m=bandrange(A)[end]
    Q=Vector[eye(n,m)[:,k] for k=1:m]
    
    
    ret=Q[1]
     
    k=0
    
    while norm(M.data.data[k+1,:])>eps()
        k+=1
        
        for j=1:m-1
            M,a,b=givensreduceab!(M,k,k+j,k,n)   
            Q[1],Q[1+j]=a*Q[1]+b*Q[1+j],-b*Q[1]+a*Q[1+j]
        end
        
        ret=Q[1]
            
        for j=1:m-1
            Q[j]=Q[j+1]
        end
        
        M,a,b=givensreduceab!(M,k,k+m,k,n)   
        
        Q[m] = -b*ret+a*eye(n)[:,k+m]
    end
    
    
    Fun(Q[1][1:k],domain(A))
end

