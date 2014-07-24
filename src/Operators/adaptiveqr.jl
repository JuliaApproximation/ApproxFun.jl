

export adaptiveqr!


##TODO: decide adaptive resize
function givensreduce!{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},v::Vector,k1::Integer,k2::Integer,j1::Integer)
    a=datagetindex(B,k1,j1)
    b=datagetindex(B,k2,j1)
    
    if b == 0.
        return B;
    end    
    

    
    sq=sqrt(abs2(a) + abs2(b))    
    a=a/sq;b=b/sq
    
    v[k1],v[k2] = a*v[k1] + b*v[k2],-b*v[k1] + a*v[k2]    
    
    
    #TODO: Assuming that left rows are already zero
    
    ir1=indexrange(B,k1)::Range1{Int64}
    ir2=indexrange(B,k2)::Range1{Int64}    
    
    for j = j1:ir1[end]
        B1 = datagetindex(B,k1,j)
        B2 = datagetindex(B,k2,j)
        
        fastsetindex!(B,a*B1 + b*B2,k1,j)
        fastsetindex!(B,-b*B1 + a*B2,k2,j)        
    end
    
    for j=ir1[end]+1:ir2[end]
        B1 = fillgetindex(B,k1,j)
        B2 = datagetindex(B,k2,j)
        
        fastsetindex!(B,a*B2 - b*B1,k2,j)
    end
    
    for j=1:B.numbcs
        B1 = getfilldata(B,k1,j)
        B2 = getfilldata(B,k2,j)
    
        setfilldata!(B, a*B1 + b*B2,k1,j)
        setfilldata!(B,-b*B1 + a*B2,k2,j)    
    end
    

    B
end

function givensreduce!(B::MutableAlmostBandedOperator,v::Vector,k1::Range1,j1::Integer)
    if length(k1)>1
        for k=k1[2]:k1[end]
            givensreduce!(B,v,k1[1],k,j1)
        end
    else
        B 
    end
end

givensreduce!(B::MutableAlmostBandedOperator,v::Vector,j::Integer)=givensreduce!(B,v,j:(j-bandrange(B)[1]),j)


function backsubstitution!{T<:Number}(B::MutableAlmostBandedOperator,u::Vector{T})
    n=length(u)
    b=bandrange(B)[end]::Int
    nbc = B.numbcs
    
    
    for k=n:-1:max(1,n-b)
        for j=k+1:n
            u[k]-=B[k,j]*u[j]
        end
          
        u[k] /= B[k,k]
    end
    
    pk = zeros(T,nbc)
    for k=n-b-1:-1:1
        for j=1:nbc
            pk[j] += u[k+b+1]*B.bc[j][k+b+1]
        end
        
        for j=k+1:k+b
            u[k]-=B[k,j]*u[j]
        end
        
        for j=1:nbc
            u[k] -= getfilldata(B,k,j)*pk[j]
        end
          
        u[k] /= B[k,k]
    end
  
  u
end


adaptiveqr(M,b)=adaptiveqr(M,b,eps())
adaptiveqr(M,b,tol)=adaptiveqr(M,b,tol,Inf)
adaptiveqr!(B,v,tol)=adaptiveqr!(B,v,tol,Inf)


promote_rule2{T<:Complex}(::Type{T},::Type{Complex{Float64}})=Complex{Float64}
promote_rule2{T<:Real}(::Type{T},::Type{Complex{Float64}})=Complex{Float64}
promote_rule2{T<:Real}(::Type{Complex{Float64}},::Type{T})=Complex{Float64}
promote_rule2{T<:Real}(::Type{T},::Type{Float64})=Float64
promote_rule2{T<:Integer}(::Type{Float64},::Type{T})=Float64



convertvec{T<:Number,V<:Number}(::BandedOperator{T},v::Vector{V})=convert(Vector{promote_rule2(T,V)},v)


function slnorm(u::Vector,r::Range)
    ret = 0.
    for j=r
        ret=max(abs(u[j]),ret)
    end
    ret
end

adaptiveqr{T<:Operator,V<:Number}(B::Vector{T},v::Vector{V},tol::Float64,N) = adaptiveqr!(MutableAlmostBandedOperator(B),convertvec(B[end],v),tol,N)  #May need to copy v in the future
function adaptiveqr!{V<:Number}(B::MutableAlmostBandedOperator,v::Vector{V},tol::Float64,N)  
    b=-bandrange(B)[1]
    m=100+b

    u=[v,zeros(V,m)]::Vector{V}
    
    l = length(v) + m  
    resizedata!(B,l)
    
    
    j=1
    ##TODO: we can allow early convergence
    while j <= N && (slnorm(u,j:j+b-1) > tol  || j <= length(v))
        if j + b == l
            u = [u,zeros(V,l)]::Vector{V}
            l *= 2
            resizedata!(B,l)
        end
        
        
        givensreduce!(B,u,j)
        j+=1
    end
  
    if j >= N
        warn("Maximum number of iterations " * string(N) * " reached")
    end
      
    backsubstitution!(B,u[1:max(j-1,length(v))])
end


