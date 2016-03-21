

export adaptiveqr!


# ca gives conj(a), cb gives conj(b), mb gives -b
# This is written to support matrix value and a and b

# applygivens! considers three cases:
#  The case where two entries in the banded part are manipulated
#  The case where one enty in the band part and one in the filled in part are manipulated
#  The case where the fill is manipulated
#

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


# this applygivens! is used to update fill.data
function applygivens!(ca,cb,mb,a,B::Matrix,k1::Integer,k2::Integer)
    for j = 1:size(B,2)
        @inbounds B1 = B[k1,j]
        @inbounds B2 = B[k2,j]

        @inbounds B[k1,j],B[k2,j]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end

hypot2(a,b) = hypot(a,b)
hypot2{T<:AbstractFloat}(a::Complex{T},b::Complex{T}) = hypot(hypot(a.re,b.re),hypot(a.im,b.im))
hypot2{T<:AbstractFloat}(a::T,b::Complex{T}) = hypot(hypot(a,b.re),b.im)
hypot2{T<:AbstractFloat}(a::Complex{T},b::T) = hypot(hypot(a.re,b),a.im)

function givensmatrix{T<:AbstractFloat}(a::T,b::T)
    c,s,r = Base.LinAlg.givensAlgorithm(a,b)
    c,s,-s,c
end

function givensmatrix{S<:Number,V<:Number}(a::S,b::V)
    if abs2(b) < 100eps2(V)^2
        #Warning: This is inconsistent for the case where a is negative
        return one(a),zero(b),zero(b),one(a)
    end

    r=hypot2(a,b)
    c,s=a/r,b/r
    conj(c),conj(s),-s,c
end


function givensreduceab!{T,M,R}(B::MutableOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer)
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

function givensreduce!{T,M,R}(B::MutableOperator{T,M,R},v::Array,k1::Integer,k2::Integer,j1::Integer)
    ca,cb,mb,a=givensreduceab!(B,k1,k2,j1)

    @simd for j=1:size(v,2)
        #@inbounds
        v[k1,j],v[k2,j] = ca*v[k1,j] + cb*v[k2,j],mb*v[k1,j] + a*v[k2,j]
    end

    B
end

function givensreduce!(B::MutableOperator,v::Array,k1::Range,j1::Integer)
    if length(k1)>1
        for k=k1[2]:k1[end]
            givensreduce!(B,v,k1[1],k,j1)
        end
    end
    B
end

givensreduce!(B::MutableOperator,v::Array,j::Integer)=givensreduce!(B,v,j:(j-bandinds(B)[1]),j)


function backsubstitution!(B::MutableOperator,u::Array)
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
adaptiveqr!(B,v)=adaptiveqr!(B,v,eps())
adaptiveqr!(B,v,tol)=adaptiveqr!(B,v,tol,Inf)






convertvec{T<:Number,V<:Number,k}(::Operator{T},v::Array{V,k})=convert(Array{promote_type(T,V),k},v)
convertvec{T<:AbstractMatrix,V<:Number,k}(::BandedOperator{T},v::Array{V,k})=totree(v)
convertvec{T<:AbstractMatrix,V<:AbstractArray,k}(::BandedOperator{T},v::Array{V,k})=convert(Array{V(promote_type(eltype(T),eltype(V))),k},v)

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
adaptiveqr{T<:Operator}(B::Vector{T},v::Array,tol::Real,N) = adaptiveqr!(MutableOperator(B),convertvec(B[end],v),tol,N)  #May need to copy v in the future
function adaptiveqr!(B::MutableOperator,v::Array,tol::Real,N)
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
        warn("Maximum number of iterations $N reached without achieving tolerance $tol.  Check that the correct number of boundary conditions are specified, or change maxiteration.")
    end

    ##TODO: why max original length?
    backsubstitution!(B,isa(u,Vector)?u[1:max(j-1,length(v))]:u[1:max(j-1,size(v,1)),:])
end

# Apply Householder(w) to the part of the operator within the bands
function applyhouseholder!(w,B::BandedMatrix,kr::Range,jr::Range)
    @simd for j = jr
        v = slice(B,kr,j)
        BLAS.axpy!(-2*(v⋅w),w,v)   # v[:]= -2*(v⋅w)*w
    end
    B
end

# Apply Householder(w) to the part of the operator both in and out of the bands
function applyhouseholder!(w,F::FillMatrix,B::BandedMatrix,kr::Range,jr::Range)
    n = size(B.data,1)
    bnd = -bandinds(B)[1]
    kr1 = kr[1]
    kr2 = kr[1]
    kr3 = kr[1]+1
    kr4 = kr[end]

    for j = jr
        x2 = slice(B,kr3:kr4,j)

        #dt represents w ⋅ [F[kr1:kr2,j];x2]
        dt=zero(eltype(w))
        m=(kr2-kr1+1)
        @inbounds for k = 1:m
            dt += w[k]⋅F[k+kr1-1,j]
        end
        @inbounds for k = 1:length(x2)
            dt+= w[k+m]⋅x2[k]
        end

        # now update x2:
        for k = eachindex(x2)
            x2[k] -= 2*dt*w[k+m]
        end
        kr2 += 1
        kr3 += 1
    end
    B
end

# Apply Householder(w) to the part of the operator outside the bands
function applyhouseholder!(w,B::Matrix,kr::Range)
    k1 = first(kr)
    m = length(w)-1


    for j = 1:size(B,2)
        dt=slice(B,k1:k1+m,j)⋅w
        @inbounds for k=1:m+1
            B[k+k1-1,j] -= 2dt*w[k]
        end
    end
    B
end

function householderreducevec!{T,M,R}(B::MutableOperator{T,M,R},kr::Range,j1::Integer, w)
    bnd=bandinds(B)[end]
    w[:] = B[kr,j1]
    w[1]-= norm(w)
    w /= norm(w)
    applyhouseholder!(w,B.data,kr,j1:kr[1]+bnd)
    applyhouseholder!(w,B.fill,B.data,kr,kr[1]+bnd+1:kr[end]+bnd)
    applyhouseholder!(w,B.fill.data,kr)
end

function householderreduce!(B::MutableOperator,v::Array,kr::Range,j1::Integer,w)
    householderreducevec!(B,kr,j1,w)

    m=length(w)
    k1=first(kr)
    for j=1:size(v,2)
        dt=slice(v,kr,j)⋅w
        @simd for k=1:m
            @inbounds v[k+k1-1,j] -= 2*dt*w[k]
        end
    end
    B
end

householderreduce!(B::MutableOperator,v::Array,j::Integer,w)=householderreduce!(B,v,j:(j-bandinds(B)[1]),j,w)

function householderadaptiveqr!(B::MutableOperator,v::Array,tol::Real,N)
    b=-B.bandinds[1]
    w = zeros(b+1) # The vector in Householder matrix is only allocated once
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

        householderreduce!(B,u,j,w)
        j+=1
    end

    if j >= N
        warn("Maximum number of iterations $N reached without achieving tolerance $tol.  Check that the correct number of boundary conditions are specified, or change maxiteration.")
    end

    ##TODO: why max original length?
    backsubstitution!(B,isa(u,Vector)?u[1:max(j-1,length(v))]:u[1:max(j-1,size(v,1)),:])
end
