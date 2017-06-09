##
# LowRankFun represents f(x,y) by r column and row slices stored as Funs.
##

export LowRankFun

"""
`LowRankFun` gives an approximation to a bivariate function in low rank form.
"""

type LowRankFun{S<:Space,M<:Space,SS<:AbstractProductSpace,T<:Number} <: BivariateFun{T}
    A::Vector{VFun{S,T}}
    B::Vector{VFun{M,T}}
    space::SS

    function LowRankFun{S,M,SS,T}(A::Vector{VFun{S,T}},
                                  B::Vector{VFun{M,T}},
                                  space::SS) where {S,M,SS,T}
        @assert length(A) == length(B)
        @assert length(A) > 0
        new{S,M,SS,T}(A,B,space)
    end
end


LowRankFun{S,M,SS,T}(A::Vector{VFun{S,T}},B::Vector{VFun{M,T}},space::SS) =
    LowRankFun{S,M,SS,T}(A,B,space)
LowRankFun{S,M,T}(A::Vector{VFun{S,T}},B::Vector{VFun{M,T}}) =
    LowRankFun(A,B,space(first(A))⊗space(first(B)))
LowRankFun{S,M,T,V}(A::Vector{VFun{S,T}},B::Vector{VFun{M,V}}) =
    LowRankFun(convert(Vector{VFun{S,promote_type(T,V)}},A),
               convert(Vector{VFun{M,promote_type(T,V)}},B),
               space(first(A))⊗space(first(B)))
Base.rank(f::LowRankFun) = length(f.A)
Base.size(f::LowRankFun,k::Integer) = k==1?mapreduce(length,max,f.A):mapreduce(length,max,f.B)
Base.size(f::LowRankFun) = size(f,1),size(f,2)

## Construction via a Matrix of coefficients

function LowRankFun{S<:Space,M<:Space,T<:Number}(X::Array{T},dx::S,dy::M)
    U,Σ,V=svd(X)
    m=max(1,count(s->s>10eps(T),Σ))

    A=VFun{S,T}[Fun(dx,U[:,k].*sqrt(Σ[k])) for k=1:m]
    B=VFun{M,T}[Fun(dy,conj(V[:,k]).*sqrt(Σ[k])) for k=1:m]

    LowRankFun(A,B)
end

## Construction in a TensorSpace via a Vector of Funs

function LowRankFun(X::Vector{VFun{S,T}},d::TensorSpace{SV,DD}) where {S,T,DD<:BivariateDomain,SV}
    @assert d[1] == space(X[1])
    LowRankFun(X,d[2])
end

function LowRankFun{S,T}(X::Vector{VFun{S,T}},dy::Space)
    m=mapreduce(ncoefficients,max,X)
    M=zeros(T,m,length(X))
    for k=1:length(X)
        M[1:ncoefficients(X[k]),k]=X[k].coefficients
    end

    LowRankFun(M,space(X[1]),dy)
end


## Adaptive constructor selector

function LowRankFun(f::F,dx::Space,dy::Space;
                    method::Symbol=:standard,tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,
                    retmax::Bool=false,gridx::Integer=64,gridy::Integer=64,maxrank::Integer=100)
    if method == :standard
        F,maxabsf=standardLowRankFun(f,dx,dy;tolerance=tolerance,gridx=gridx,gridy=gridy,maxrank=maxrank)
    elseif method == :Cholesky
        @assert domain(dx) == domain(dy)
        if dx == dy
            F,maxabsf=CholeskyLowRankFun(f,dx;tolerance=tolerance,grid=max(gridx,gridy),maxrank=maxrank)
        else
            G,maxabsf=CholeskyLowRankFun(f,dx;tolerance=tolerance,grid=max(gridx,gridy),maxrank=maxrank)
            F=LowRankFun(copy(G.A),map(b->Fun(b,dy),G.B))
        end
    end
    retmax ? (F,maxabsf) : F
end

## Standard adaptive construction

function standardLowRankFun(f::F,dx::Space,dy::Space;tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,gridx::Integer=64,gridy::Integer=64,maxrank::Integer=100)
    xy = checkpoints(dx⊗dy)
    T = promote_type(eltype(f(first(xy)...)),prectype(dx),prectype(dy))

    # We start by sampling on the given grid, find the approximate maximum and create the first rank-one approximation.
    ptsx,ptsy=points(dx,gridx),points(dy,gridy)
    X = zeros(T,gridx,gridy)
    maxabsf,r=findapproxmax!(f,X,ptsx,ptsy,gridx,gridy)
    if maxabsf < eps(zero(T))/eps(T) return LowRankFun([Fun(dx,[zero(T)])],[Fun(dy,[zero(T)])]),maxabsf end
    a,b=Fun(x->f(x,r[2]),dx),Fun(y->f(r[1],y),dy)

    # If necessary, we resize the grid to be at least as large as the
    # lengths of the first row and column Funs and we recompute the values of X.
    if gridx < ncoefficients(a) || gridy < ncoefficients(b)
        gridx,gridy = max(gridx,ncoefficients(a)),max(gridy,ncoefficients(b))
        ptsx,ptsy=points(dx,gridx),points(dy,gridy)
        X = zeros(T,gridx,gridy)
        maxabsf,r=findapproxmax!(f,X,ptsx,ptsy,gridx,gridy)
        a,b=Fun(x->f(x,r[2]),dx),Fun(y->f(r[1],y),dy)
    end

    A,B=typeof(a)[],typeof(b)[]
    if tolerance == :relative
        tol = 100maxabsf*eps(T)
    elseif tolerance[1] == :absolute
        tol = 100*tolerance[2]*eps(T)
    end
    tol10 = tol/10
    Avals,Bvals = zeros(T,gridx),zeros(T,gridy)
    p₁,p₂ = plan_transform(dx,Avals),plan_transform(dy,Bvals)

    # Eat, drink, subtract rank-one, repeat.
    for k=1:maxrank
        if (norm(a.coefficients,Inf) < tol || norm(b.coefficients,Inf) < tol)
            return LowRankFun(A,B),maxabsf
        end
        A,B =[A;a/sqrt(abs(a(r[1])))],[B;b/(sqrt(abs(b(r[2])))*sign(b(r[2])))]
        r=findapproxmax!(A[k],B[k],X,ptsx,ptsy,gridx,gridy)
        Ar,Br=evaluate(A,r[1]),evaluate(B,r[2])
        for i=1:gridx
            @inbounds Avals[i] = f(ptsx[i],r[2])
        end
        for j=1:gridy
            @inbounds Bvals[j] = f(r[1],ptsy[j])
        end
        a,b = Fun(dx,p₁*Avals) - dotu(Br,A),Fun(dy,p₂*Bvals) - dotu(Ar,B)
        chop!(a,tol10),chop!(b,tol10)
    end
    warn("Maximum rank of " * string(maxrank) * " reached")
    return LowRankFun(A,B),maxabsf
end

## Adaptive Cholesky decomposition, when f is Hermitian positive (negative) definite

function CholeskyLowRankFun(f::F,dx::Space;tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,grid::Integer=64,maxrank::Integer=100)
    xy = checkpoints(dx⊗dx)
    T = promote_type(eltype(f(first(xy)...)),prectype(dx))

    # We start by sampling on the given grid, find the approximate maximum and create the first rank-one approximation.
    pts=points(dx,grid)
    X = zeros(T,grid)
    maxabsf,r=findcholeskyapproxmax!(f,X,pts,grid)
    if maxabsf < eps(zero(T))/eps(T) return LowRankFun([Fun(dx,[zero(T)])],[Fun(dx,[zero(T)])]),maxabsf end
    a=Fun(x->f(x,r),dx)

    # If necessary, we resize the grid to be at least as large as the
    # ncoefficients of the first row/column Fun and we recompute the values of X.
    if grid < ncoefficients(a)
        grid = max(grid,ncoefficients(a))
        pts=points(dx,grid)
        X = zeros(T,grid)
        maxabsf,r=findcholeskyapproxmax!(f,X,pts,grid)
        a=Fun(x->f(x,r),dx)
    end

    A,B=typeof(a)[],typeof(a)[]
    if tolerance == :relative
        tol = 100maxabsf*eps(T)
    elseif tolerance[1] == :absolute
        tol = 100*tolerance[2]*eps(T)
    end
    tol10 = tol/10
    Avals = zeros(T,grid)
    p₁ = plan_transform(dx,Avals)

    # Eat, drink, subtract rank-one, repeat.
    for k=1:maxrank
        if norm(a.coefficients,Inf) < tol return LowRankFun(A,B),maxabsf end
        A,B = [A;a/sqrt(abs(a(r)))],[B;a/(sqrt(abs(a(r)))*sign(a(r)))]
        r=findcholeskyapproxmax!(A[k],B[k],X,pts,grid)
        Br=evaluate(B,r)
        for i=1:grid
            @inbounds Avals[i] = f(pts[i],r)
        end
        a = Fun(dx,p₁*Avals) - dotu(Br,A)
        chop!(a,tol10)
    end
    warn("Maximum rank of " * string(maxrank) * " reached")
    return LowRankFun(A,B),maxabsf
end


## Construction via TensorSpaces and ProductDomains

LowRankFun(f::Function,args...;kwds...) = LowRankFun(F(f),args...;kwds...)

LowRankFun(f::F,S::TensorSpace{SV,DD,RR};kwds...) where {SV,DD<:BivariateDomain,RR} =
    LowRankFun(f,S[1],S[2];kwds...)
LowRankFun(f::F,dx::Domain,dy::Domain;kwds...) =
    LowRankFun(f,Space(dx),Space(dy);kwds...)
LowRankFun(f::F,d::ProductDomain;kwds...) =
    LowRankFun(f,d.domains...;kwds...)

LowRankFun(f::F,d1::Vector,d2::Vector;kwds...) =
    LowRankFun(f,convert(Domain,d1),convert(Domain,d2);kwds...)
LowRankFun(f::F;kwds...) = LowRankFun(f,Interval(),Interval();kwds...)

## Construction from values

LowRankFun{T<:Number}(A::Array{T}) = LowRankFun(A,Interval{T}(),Interval{T}())
LowRankFun(c::Number,etc...) = LowRankFun((x,y)->c,etc...)

## Construction from other LowRankFuns

LowRankFun(f::LowRankFun,d1::IntervalDomain,d2::IntervalDomain) =
    LowRankFun(map(g->Fun(d1,g.coefficients),f.A),
               map(g->Fun(d2,g.coefficients),f.B))
LowRankFun(f::LowRankFun) = LowRankFun(f,Interval(),Interval())



## Utilities

function findapproxmax!(f::F,X::AbstractMatrix,ptsx::AbstractVector,ptsy::Vector,gridx,gridy)
    for j=1:gridy
        ptsyj = ptsy[j]
        @simd for k=1:gridx
            @inbounds X[k,j]+=f(ptsx[k],ptsyj)
        end
    end
    maxabsf,impt = findmaxabs(X)
    imptple = ind2sub((gridx,gridy),impt)
    maxabsf,[ptsx[imptple[1]],ptsy[imptple[2]]]
end

function findapproxmax!(A::Fun,B::Fun,X::AbstractMatrix,ptsx::Vector,ptsy::Vector,gridx,gridy)
    Ax,By = A.(ptsx),B.(ptsy)
    subtractrankone!(Ax,By,X,gridx,gridy)
    maxabsf,impt = findmaxabs(X)
    imptple = ind2sub((gridx,gridy),impt)
    [ptsx[imptple[1]],ptsy[imptple[2]]]
end

function findcholeskyapproxmax!(f::F,X::AbstractVector,pts::Vector,grid)
    @simd for k=1:grid
        @inbounds X[k]+=f(pts[k],pts[k])
    end
    maxabsf,impt = findmaxabs(X)
    maxabsf,pts[impt]
end

function findcholeskyapproxmax!(A::Fun,B::Fun,X::Vector,pts::Vector,grid)
    Ax,By = A.(pts),B.(pts)
    subtractrankone!(Ax,By,X,grid)
    maxabsf,impt = findmaxabs(X)
    pts[impt]
end

function subtractrankone!(A::AbstractVector,B::AbstractVector,X::AbstractMatrix,gridx::Int,gridy::Int)
    for j=1:gridy
        @inbounds Bj = B[j]
        @simd for k=1:gridx
            @inbounds X[k,j] -= A[k]*Bj
        end
    end
end

function subtractrankone!(A::AbstractVector,B::AbstractVector,X::AbstractVector,grid::Int)
    @simd for k=1:grid
        @inbounds X[k] -= A[k]*B[k]
    end
end

## TODO: in Julia base?
function findmaxabs(a)
    if isempty(a)
        throw(ArgumentError("collection must be non-empty"))
    end
    m = abs(a[1])
    mi = 1
    for i in eachindex(a)
        ai = abs(a[i])
        if ai > m || m!=m
            m = ai
            mi = i
        end
    end
    return (m, mi)
end

(f::LowRankFun)(x,y)=evaluate(f,x,y)

domain(f::LowRankFun,k::Integer)=k==1? domain(first(f.A)) : domain(first(f.B))
space(f::LowRankFun,k::Integer)=k==1? space(first(f.A)) : space(first(f.B))
space(f::LowRankFun)=f.space

Base.transpose{S,M,SS,T}(f::LowRankFun{S,M,SS,T})=LowRankFun(f.B,f.A,transpose(space(f)))

function values(f::LowRankFun)
    xm=mapreduce(ncoefficients,max,f.A)
    ym=mapreduce(ncoefficients,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=values(pad(f.A[k],xm))*values(pad(f.B[k],ym)).'
    end
    ret
end

#TODO: this is inconsistent with 1D where it does canonical
function coefficients(f::LowRankFun)
    xm=mapreduce(ncoefficients,max,f.A)
    ym=mapreduce(ncoefficients,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(f.A[k].coefficients,xm)*pad(f.B[k].coefficients,ym).'
    end
    ret
end

function coefficients(f::LowRankFun,n::Space,m::Space)
    xm=mapreduce(ncoefficients,max,f.A)
    ym=mapreduce(ncoefficients,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(coefficients(f.A[k],n),xm)*pad(coefficients(f.B[k],m),ym).'
    end
    ret
end

function vecpoints(f::LowRankFun,k::Integer)
    if k==1
        xm=mapreduce(ncoefficients,max,f.A)
        points(space(first(f.A)),xm)
    else
        ym=mapreduce(ncoefficients,max,f.B)
        points(space(first(f.B)),ym)
    end
end


evaluate{T<:Fun,M<:Fun}(A::Vector{T},B::Vector{M},x,y)=dotu(evaluate(A,x),evaluate(B,y))
evaluate{T<:Fun,M<:Fun}(A::Vector{T},B::Vector{M},x::AbstractVector,y::AbstractVector)=evaluate.(A.',x)*evaluate.(B,y.')

evaluate(f::LowRankFun,x,y)=evaluate(f.A,f.B,x,y)
evaluate(f::LowRankFun,::Colon,::Colon)=f
evaluate(f::LowRankFun,x::Number,::Colon)=dotu(f.B,evaluate(f.A,x))

function evaluate(f::LowRankFun,::Colon,y::Number)
    m = maximum(map(ncoefficients,f.A))
    r=rank(f)
    ret = zeros(m)

    for k=1:r
        for j=1:ncoefficients(f.A[k])
            @inbounds ret[j] += f.A[k].coefficients[j]*f.B[k](y)
        end
    end

    Fun(first(f.A).space,ret)
end


## Truncate
#TODO: should reduce rank if needed
Base.chop(f::LowRankFun,tol)=LowRankFun(map(g->chop(g,tol),f.A),map(g->chop(g,tol),f.B),f.space)
function pad(f::LowRankFun,m::Integer,n::Integer)
    A,B = deepcopy(f.A),deepcopy(f.B)

    for k=1:rank(f)
        pad!(A[k],m);pad!(B[k],n)
    end

    LowRankFun(A,B,f.space)
end


## Algebra

for op = (:*,:/)
    @eval ($op){T<:Fun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::LowRankFun,c::Number) = LowRankFun(($op)(f.A,c),f.B)
    @eval ($op)(c::Number,f::LowRankFun) = LowRankFun(($op)(c,f.A),f.B)
end

# Let K be a LowRankFun and f be a Fun.
# op(f,K) acts as operating in the x variable, and
# op(K,f) acts as operating in the y variable.

for op = (:*,:/)
    @eval ($op)(f::Fun,K::LowRankFun) = LowRankFun(($op)(f,K.A),K.B)
    @eval ($op)(K::LowRankFun,f::Fun) = LowRankFun(K.A,($op)(K.B,f))
end

+(f::LowRankFun,g::LowRankFun) = LowRankFun([f.A;g.A],[f.B;g.B])
-(f::LowRankFun) = LowRankFun(-f.A,f.B)
-(f::LowRankFun,g::LowRankFun) = f+(-g)

## QR factorization of a LowRankFun

function Base.qr(f::LowRankFun)
    sp,r = space(f),rank(f)
    Q,R = qr(coefficients(f.A))
    BR = coefficients(f.B)*R.'
    LowRankFun(map(i->Fun(sp[1],Q[:,i]),1:r),map(i->Fun(sp[2],BR[:,i]),1:r),sp)
end

## Special functions

Base.real(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(real,u.B),-map(imag,u.B)])
Base.imag(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(imag,u.B),map(real,u.B)])


## Calculus


Base.sum(g::LowRankFun)=dotu(map(sum,g.A),map(sum,g.B))
Base.sum(g::LowRankFun,n::Integer)=(n==1)?dotu(g.B,map(sum,g.A)):dotu(g.A,map(sum,g.B))
Base.cumsum(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(cumsum,g.A),copy(g.B)):LowRankFun(copy(g.A),map(cumsum,g.B))
differentiate(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(differentiate,g.A),copy(g.B)):LowRankFun(copy(g.A),map(differentiate,g.B))
integrate(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(integrate,g.A),copy(g.B)):LowRankFun(copy(g.A),map(integrate,g.B))
