

export adaptiveqr!

function applygivens!(a,b,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    ca,cb=conj(a),conj(b)
    @simd for j = jr
        @inbounds B1 = B.data[j-k1+B.l+1,k1]    #B[k1,j]
        @inbounds B2 = B.data[j-k2+B.l+1,k2]    #B[k2,j]       
        
        @inbounds B.data[j-k1+B.l+1,k1],B.data[j-k2+B.l+1,k2]= ca*B1 + cb*B2,-b*B1 + a*B2
    end   
    
    B
end

function applygivens!(a,b,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        B1 = unsafe_getindex(F,k1,j)
        @inbounds B2 = B.data[j-k2+B.l+1,k2]   #B[k2,j]     
        
        @inbounds B.data[j-k2+B.l+1,k2]=a*B2 - b*B1
    end   
    
    B
end

function applygivens!(a,b,B::Matrix,k1::Integer,k2::Integer)
    ca,cb=conj(a),conj(b)
    for j = 1:size(B,2)
        @inbounds B1 = B[k1,j]
        @inbounds B2 = B[k2,j]
        
        @inbounds B[k1,j],B[k2,j]= ca*B1 + cb*B2,-b*B1 + a*B2
    end   
    
    B
end


function givensreduceab!{T<:Number,M,R}(B::AlmostBandedOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer)
    bnd=B.bandinds
    A=B.data
    
    @inbounds a=A.data[j1-k1+A.l+1,k1]  #A[k1,j1]
    @inbounds b=A.data[j1-k2+A.l+1,k2]  #A[k2,j1] 
    
    if b == 0
        return one(T),zero(T)
    end    
    
    sq=sqrt(abs2(a) + abs2(b))    
    a=a/sq;b=b/sq
    
    
    #Assuming that left rows are already zero    
    applygivens!(a,b,B.data,k1,k2,j1:k1+B.bandinds[end])
    applygivens!(a,b,B.fill,B.data,k1,k2,k1+B.bandinds[end]+1:k2+B.bandinds[end])    
    applygivens!(a,b,B.fill.data,k1,k2)


    a::T,b::T
end

function givensreduce!{T<:Number,M,R}(B::AlmostBandedOperator{T,M,R},v::Array,k1::Integer,k2::Integer,j1::Integer)
    a,b=givensreduceab!(B,k1,k2,j1)
    
    if b != 0.0
        ca=conj(a)
        cb=conj(b)
    
        @simd for j=1:size(v,2)
            #@inbounds 
            v[k1,j],v[k2,j] = ca*v[k1,j] + cb*v[k2,j],-b*v[k1,j] + a*v[k2,j]    
        end
    end        

    B
end

function givensreduce!(B::AlmostBandedOperator,v::Array,k1::UnitRange,j1::Integer)
    if length(k1)>1
        for k=k1[2]:k1[end]
            givensreduce!(B,v,k1[1],k,j1)
        end
    end
    B 
end

givensreduce!(B::AlmostBandedOperator,v::Array,j::Integer)=givensreduce!(B,v,j:(j-bandinds(B)[1]),j)


function backsubstitution!{T<:Number}(B::AlmostBandedOperator,u::Array{T})
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.fill.numbcs
    A=B.data
    
    pk = zeros(T,nbc)
    
    for c=1:size(u,2)
        fill!(pk,zero(T))
    
        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            @simd for j=k+1:n             
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end
              
            @inbounds u[k,c] /= A.data[A.l+1,k]
        end
        
       #filled rows
        for k=n-b-1:-1:1
            @simd for j=1:nbc       
                @inbounds pk[j] += u[k+b+1,c]*B.fill.bc[j].data[k+b+1]
            end
            
            @simd for j=k+1:k+b          
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end
            
            @simd for j=1:nbc        
                @inbounds u[k,c] -= B.fill.data[k,j]*pk[j]
            end
              
            @inbounds u[k,c] /= A.data[A.l+1,k]
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
   for k=r
        @simd for j=1:size(u,2)
            #@inbounds 
            ret=max(abs(u[k,j]),ret)
        end
    end
    ret
end

function slnorm(u::BandedMatrix,r::Range)
    ret = 0.0
   for k=r
        @simd for j=1:size(u.data,1)
            #@inbounds 
            ret=max(abs(u.data[j,k]),ret)
        end
    end
    ret
end

adaptiveqr{V<:Number}(B::Operator,v::Array{V},tol::Float64,N) = adaptiveqr([B],v,tol,N)  #May need to copy v in the future
adaptiveqr{T<:Operator,V<:Number}(B::Vector{T},v::Array{V},tol::Float64,N) = adaptiveqr!(AlmostBandedOperator(B),convertvec(B[end],v),tol,N)  #May need to copy v in the future
function adaptiveqr!{V<:Number}(B::AlmostBandedOperator,v::Array{V},tol::Float64,N)  
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


