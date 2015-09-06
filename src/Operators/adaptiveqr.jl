

export adaptiveqr!

function applygivens!(ca,cb,mb,a,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    @simd for j = jr
        @inbounds B1 = B.data[j-k1+B.l+1,k1]    #B[k1,j]
        @inbounds B2 = B.data[j-k2+B.l+1,k2]    #B[k2,j]

        @inbounds B.data[j-k1+B.l+1,k1],B.data[j-k2+B.l+1,k2]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end

function applygivens!(ca,cb,mb,a,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        B1 = unsafe_getindex(F,k1,j)
        @inbounds B2 = B.data[j-k2+B.l+1,k2]   #B[k2,j]

        @inbounds B.data[j-k2+B.l+1,k2]=mb*B1 + a*B2
    end

    B
end

function applygivens!(ca,cb,mb,a,B::Matrix,k1::Integer,k2::Integer)
    for j = 1:size(B,2)
        @inbounds B1 = B[k1,j]
        @inbounds B2 = B[k2,j]

        @inbounds B[k1,j],B[k2,j]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end


function givensmatrix(a::Number,b::Number)
    if abs(b) < 10eps()
        #Warning: This is inconsistent for the case where a is negative
        return one(a),zero(b),zero(b),one(a)
    end

    sq=sqrt(abs2(a) + abs2(b))
    a,b=a/sq,b/sq
    conj(a),conj(b),-b,a
end


function givensreduceab!{T,M,R}(B::AlmostBandedOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer)
    bnd=B.bandinds
    A=B.data

    @inbounds a=A.data[j1-k1+A.l+1,k1]  #A[k1,j1]
    @inbounds b=A.data[j1-k2+A.l+1,k2]  #A[k2,j1]

    ca,cb,mb,a=givensmatrix(a,b)


    #Assuming that left rows are already zero
    applygivens!(ca,cb,mb,a,B.data,k1,k2,j1:k1+B.bandinds[end])
    applygivens!(ca,cb,mb,a,B.fill,B.data,k1,k2,k1+B.bandinds[end]+1:k2+B.bandinds[end])
    applygivens!(ca,cb,mb,a,B.fill.data,k1,k2)


    ca,cb,mb,a
end

function givensreduce!{T,M,R}(B::AlmostBandedOperator{T,M,R},v::Array,k1::Integer,k2::Integer,j1::Integer)
    ca,cb,mb,a=givensreduceab!(B,k1,k2,j1)

    @simd for j=1:size(v,2)
        #@inbounds
        v[k1,j],v[k2,j] = ca*v[k1,j] + cb*v[k2,j],mb*v[k1,j] + a*v[k2,j]
    end

    B
end

function givensreduce!(B::AlmostBandedOperator,v::Array,k1::Range,j1::Integer)
    if length(k1)>1
        for k=k1[2]:k1[end]
            givensreduce!(B,v,k1[1],k,j1)
        end
    end
    B
end

givensreduce!(B::AlmostBandedOperator,v::Array,j::Integer)=givensreduce!(B,v,j:(j-bandinds(B)[1]),j)


function backsubstitution!(B::AlmostBandedOperator,u::Array)
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.fill.numbcs
    A=B.data
    T=eltype(u)
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






convertvec{T<:Number,V<:Number,k}(::BandedOperator{T},v::Array{V,k})=convert(Array{promote_type(T,V),k},v)
convertvec{T<:AbstractMatrix,V<:Number,k}(::BandedOperator{T},v::Array{V,k})=totree(v)
convertvec{T<:AbstractMatrix,V<:AbstractArray,k}(::BandedOperator{T},v::Array{V,k})=convert(Array{V(promote_type(eltype(T),eltype(V))),k},v)

function slnorm{T}(u::Array{T},r::Range)
    ret = zero(real(T))
    for k=r
        @simd for j=1:size(u,2)
            #@inbounds
            ret=max(norm(u[k,j]),ret)
        end
    end
    ret
end

function slnorm{T}(u::BandedMatrix{T},r::Range)
    ret = zero(real(T))
    for k=r
        @simd for j=1:size(u.data,1)
            #@inbounds
            ret=max(abs(u.data[j,k]),ret)
        end
    end
    ret
end

adaptiveqr(B::Operator,v::Array,tol::Real,N) = adaptiveqr([B],v,tol,N)  #May need to copy v in the future
adaptiveqr{T<:Operator}(B::Vector{T},v::Array,tol::Real,N) = adaptiveqr!(AlmostBandedOperator(B),convertvec(B[end],v),tol,N)  #May need to copy v in the future
function adaptiveqr!(B::AlmostBandedOperator,v::Array,tol::Real,N)
    b=-B.bandinds[1]
    m=3+b

    l = size(v,1) + m

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
        warn("Maximum number of iterations " * string(N) * " reached.  Check that the correct number of boundary conditions are specified, or change maxiteration.")
    end

    ##TODO: why max original length?
    backsubstitution!(B,isa(u,Vector)?u[1:max(j-1,length(v))]:u[1:max(j-1,size(v,1)),:])
end
