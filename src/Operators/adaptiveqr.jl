

export adaptiveqr!


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
    
    for j=1:B.numbcs
        B1 = getfilldata(B,k1,j)
        B2 = getfilldata(B,k2,j)
    
        setfilldata!(B, a*B1 + b*B2,k1,j)
        setfilldata!(B,-b*B1 + a*B2,k2,j)    
    end
    

    a::T,b::T
end

function givensreduce!{T<:Number,M,R}(B::MutableAlmostBandedOperator{T,M,R},v::Array,k1::Integer,k2::Integer,j1::Integer)
    a,b=givensreduceab!(B,k1,k2,j1)
    
    if b != 0.0
        for j=1:size(v,2)
            v[k1,j],v[k2,j] = a*v[k1,j] + b*v[k2,j],-b*v[k1,j] + a*v[k2,j]    
        end
    end        

    B
end

function givensreduce!(B::MutableAlmostBandedOperator,v::Array,k1::Range1,j1::Integer)
    if length(k1)>1
        for k=k1[2]:k1[end]
            givensreduce!(B,v,k1[1],k,j1)
        end
    else
        B 
    end
end

givensreduce!(B::MutableAlmostBandedOperator,v::Array,j::Integer)=givensreduce!(B,v,j:(j-bandrange(B)[1]),j)


function backsubstitution!{T<:Number}(B::MutableAlmostBandedOperator,u::Array{T})
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.numbcs
    
     pk = zeros(T,nbc)
    
    for c=1:size(u,2)
        fill!(pk,zero(T))
    
        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            for j=k+1:n
                u[k,c]-=B[k,j]*u[j,c]
            end
              
            u[k,c] /= B[k,k]
        end
        
       #filled rows
        for k=n-b-1:-1:1
            for j=1:nbc
                pk[j] += u[k+b+1,c]*B.bc[j][k+b+1]
            end
            
            for j=k+1:k+b
                u[k,c]-=B[k,j]*u[j,c]
            end
            
            for j=1:nbc
                u[k,c] -= getfilldata(B,k,j)*pk[j]
            end
              
            u[k,c] /= B[k,k]
        end
    end
    u
end


adaptiveqr(M,b)=adaptiveqr(M,b,eps())
adaptiveqr(M,b,tol)=adaptiveqr(M,b,tol,Inf)
adaptiveqr!(B,v,tol)=adaptiveqr!(B,v,tol,Inf)






convertvec{T<:Number,V<:Number}(::BandedOperator{T},v::Vector{V})=convert(Vector{promote_type(T,V)},v)
convertvec{T<:Number,V<:Number}(::BandedOperator{T},v::Array{V,2})=convert(Array{promote_type(T,V),2},v)

function slnorm(u::Array,r::Range)
    ret = 0.0
    for k=r,j=1:size(u,2)
        ret=max(abs(u[k,j]),ret)
    end
    ret
end
adaptiveqr{V<:Number}(B::Operator,v::Array{V},tol::Float64,N) = adaptiveqr([B],v,tol,N)  #May need to copy v in the future
adaptiveqr{T<:Operator,V<:Number}(B::Vector{T},v::Array{V},tol::Float64,N) = adaptiveqr!(MutableAlmostBandedOperator(B),convertvec(B[end],v),tol,N)  #May need to copy v in the future
function adaptiveqr!{V<:Number}(B::MutableAlmostBandedOperator,v::Array{V},tol::Float64,N)  
    b=-B.bandinds[1]
    m=100+b
    
    l = length(v) + m  
    
    u=pad(v,l,size(v,2))    
    resizedata!(B,l)
    
    
    j=1
    ##TODO: we can allow early convergence
    while j <= N && (slnorm(u,j:j+b-1) > tol  || j <= size(v,1))
        if j + b == l
            l *= 2
            u = pad(u,l,size(u,2))            
            resizedata!(B,l)
        end
        
        
        givensreduce!(B,u,j)
        j+=1
    end
  
    if j >= N
        warn("Maximum number of iterations " * string(N) * " reached")
    end
      
    ##TODO: why max original length?
    backsubstitution!(B,isa(u,Vector)?u[1:max(j-1,length(v))]:u[1:max(j-1,size(v,1)),:])
end


