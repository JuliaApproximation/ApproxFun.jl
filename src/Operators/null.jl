


#assume dimension of kernel is range space
Base.null(A::BandedOperator)=null(A,order(rangespace(A)))


function applygivens!(a,b,Q1,Qk)
    @simd for j=1:length(Q1)
        @inbounds Q1[j],Qk[j]=a*Q1[j]+b*Qk[j],-b*Q1[j]+a*Qk[j]
    end
end



#d is number of elements in the kernel

function Base.null{T<:Number}(A::BandedOperator{T},d,maxit=100000)
    M=MutableOperator(A')
    m=bandinds(A)[end]
    n=m+100
    resizedata!(M,n)
    Q=Vector{T}[zeros(T,n) for j=1:m]
    for j=1:m
        Q[j][j]=one(T)
    end

    Q1=Q[1]

    k=0

    while slnorm(M.data,k+1:k+d)>eps(T)  && k <= maxit
        k+=1

        if k+m+d >= n
            n=2n
            resizedata!(M,n)

            for j=1:m
                Q[j]=[Q[j];zeros(T,n)]
            end
        end

        Q1=Q[1]

        for j=1:m-1
            a,b=givensreduceab!(M,k,k+j,k)
            applygivens!(a,b,Q1,Q[1+j] )
        end


        for j=1:m-1
            @inbounds Q[j]=Q[j+1]
        end

        a,b=givensreduceab!(M,k,k+m,k)

        Q[m] = Q1
        Q[m] *= -b
        Q[m][k+m] += a
    end


    ds=domainspace(A)
    Fun{typeof(ds),T}[Fun(Q[j][1:k],ds) for j=1:d]
end


Base.null{T<:BandedOperator}(A::Array{T,2},d)=null(interlace(A),d)

function Base.null{T<:BandedOperator}(A::Array{T,2})
    # We assume ultraspherical space gives dimension of kernel
    d = 0
    A=promotespaces(A)
    for k=1:size(A,1)
        d += order(rangespace(A[k,1]))
    end

    null(A,d)
end


