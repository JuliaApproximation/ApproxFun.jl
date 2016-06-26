

export adaptiveqr!




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


function givensreduce!{T,M,R}(B::MutableOperator{T,M,R},v::Array,k1::Integer,k2::Integer,j1::Integer)
    ca,cb,mb,a=givensreduceab!(B,k1,k2,j1)

    @simd for j=1:size(v,2)
        @inbounds v[k1,j],v[k2,j] = ca*v[k1,j] + cb*v[k2,j],mb*v[k1,j] + a*v[k2,j]
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



adaptiveqr(M,b)=adaptiveqr(M,b,eps())
adaptiveqr(M,b,tol)=adaptiveqr(M,b,tol,Inf)
adaptiveqr!(B,v)=adaptiveqr!(B,v,eps())
adaptiveqr!(B,v,tol)=adaptiveqr!(B,v,tol,Inf)






convertvec{T<:Number,V<:Number,k}(::Operator{T},v::Array{V,k})=convert(Array{promote_type(T,V),k},v)
convertvec{T<:AbstractMatrix,V<:Number,k}(::BandedOperator{T},v::Array{V,k})=totree(v)
convertvec{T<:AbstractMatrix,V<:AbstractArray,k}(::BandedOperator{T},v::Array{V,k})=convert(Array{V(promote_type(eltype(T),eltype(V))),k},v)

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
    while j <= N && (slnorm(u,j:j+b-1,:) > tol  || j <= size(v,1))
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

function applyhouseholder!{T}(w::Vector{T},B::AbstractArray{T},k1::Int64,kr::Range,j::Int64)
    v = slice(B,k1+kr,j)
    BLAS.axpy!(-2*(v⋅w),w,v)   # v[:]= -2*(v⋅w)*w
    #dt = zero(T)
    #@inbounds @simd for k in kr dt+=B[k+k1,j]⋅w[k] end
    #dt *= -2
    #@inbounds @simd for k in kr B[k+k1,j]+=dt*w[k] end
    B
end

function applyhouseholder!{T}(w::Vector{T},B::AbstractArray{T},kr::Range,k2::Int64,j::Int64)
    v = slice(B,kr,j)
    BLAS.axpy!(-2*(v⋅w),w,v)   # v[:]= -2*(v⋅w)*w
    #dt = zero(T)
    #@inbounds @simd for k in kr dt+=B[k,j]⋅w[k+k2] end
    #dt *= -2
    #@inbounds @simd for k in kr B[k,j]+=dt*w[k+k2] end
    B
end

# Apply Householder(w) to the part of the operator within the bands
function applyhouseholder!{T}(w::Vector{T},B::BandedMatrix,kr::Range,jr::Range)
    for j = jr applyhouseholder!(w,B,kr,1-first(kr),j) end
    B
end


# Apply Householder(w) to the part of the operator both in and out of the bands
function applyhouseholder!{T}(w::Vector{T},F::FillMatrix,B::BandedMatrix,kr::Range,jr::Range)
    k1 = first(kr)-1
    k2 = last(kr)-first(kr)
    m=1

    for j = jr
        #dt represents w ⋅ [F[kr1:kr2,j];x2]
        dt=zero(T)
        @inbounds @simd for k = 1:m dt += w[k]⋅unsafe_getindex(F,k+k1,j) end
        @inbounds @simd for k = 1+m:k2 dt+=w[k]⋅B[k+k1,j] end
        dt *= -2
        @inbounds @simd for k = 1+m:k2 B[k+k1,j] += dt*w[k] end
        m+=1
    end
    B
end


# Apply Householder(w) to the part of the operator outside the bands
function applyhouseholder!{T}(w::Vector{T},B::Matrix,kr::Range)
    for j = 1:size(B,2) applyhouseholder!(w,B,first(kr)-1,1:length(w),j) end
    B
end

function householderreducevec!{T,M,R}(w::Vector{T}, B::MutableOperator{T,M,R},kr::Range,j1::Integer)
    bnd=bandinds(B)[end]
    copy!(w,B.data[kr,j1])
    w[1]-= norm(w)
    # TODO: should be normalize!(w) below, since this is compatible with any type T.
    # But this is only supported in Julia v0.5 https://github.com/JuliaLang/julia/pull/13681
    scal!(inv(norm(w)),w) #w/=norm(w)
    applyhouseholder!(w,B.data,kr,j1:kr[1]+bnd)
    applyhouseholder!(w,B.fill,B.data,kr,kr[1]+bnd+1:kr[end]+bnd)
    applyhouseholder!(w,B.fill.data,kr)
end



function householderreduce!{T}(w::Vector{T},B::MutableOperator,v::Array,kr::Range,j1::Integer)
    householderreducevec!(w,B,kr,j1)
    for j=1:size(v,2) applyhouseholder!(w,v,first(kr)-1,1:length(w),j) end
    B
end

householderreduce!(w::Vector,B::MutableOperator,v::Array,j::Integer)=householderreduce!(w,B,v,j:(j-bandinds(B)[1]),j)

function householderadaptiveqr!{T,M,R}(B::MutableOperator{T,M,R},v::Array,tol::Real,N)
    b=-B.bandinds[1]
    w = zeros(T,b+1) # The vector in Householder matrix is only allocated once
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

        householderreduce!(w,B,u,j)
        j+=1
    end

    if j >= N
        warn("Maximum number of iterations $N reached without achieving tolerance $tol.  Check that the correct number of boundary conditions are specified, or change maxiteration.")
    end

    ##TODO: why max original length?
    backsubstitution!(B,isa(u,Vector)?u[1:max(j-1,length(v))]:u[1:max(j-1,size(v,1)),:])
end
