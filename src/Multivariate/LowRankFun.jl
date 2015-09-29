##
# LowRankFun represents f(x,y) by r column and row slices stored as Funs.
##

export LowRankFun

"""
`LowRankFun` gives an approximation to a bivariate function in low rank form.
"""
type LowRankFun{S<:Space,M<:Space,SS<:AbstractProductSpace,T<:Number} <: BivariateFun{T}
    A::Vector{Fun{S,T}}
    B::Vector{Fun{M,T}}
    space::SS

    function LowRankFun(A::Vector{Fun{S,T}},B::Vector{Fun{M,T}},space::SS)
        @assert length(A) == length(B)
        @assert length(A) > 0
        new(A,B,space)
    end
end

LowRankFun{S,M,SS,T}(A::Vector{Fun{S,T}},B::Vector{Fun{M,T}},space::SS)=LowRankFun{S,M,SS,T}(A,B,space)
LowRankFun{S,M,T}(A::Vector{Fun{S,T}},B::Vector{Fun{M,T}})=LowRankFun(A,B,space(first(A))⊗space(first(B)))
LowRankFun{S,M,T,V}(A::Vector{Fun{S,T}},B::Vector{Fun{M,V}})=LowRankFun(convert(Vector{Fun{S,promote_type(T,V)}},A),convert(Vector{Fun{M,promote_type(T,V)}},B),space(first(A))⊗space(first(B)))

Base.rank(f::LowRankFun)=length(f.A)
Base.size(f::LowRankFun,k::Integer)=k==1?mapreduce(length,max,f.A):mapreduce(length,max,f.B)
Base.size(f::LowRankFun)=size(f,1),size(f,2)

## Construction via a Matrix of coefficients

function LowRankFun{S<:Space,M<:Space,T<:Number}(X::Array{T},dx::S,dy::M)
    U,Σ,V=svd(X)
    m=max(1,count(s->s>10eps(T),Σ))

    A=Fun{S,T}[Fun(U[:,k].*sqrt(Σ[k]),dx) for k=1:m]
    B=Fun{M,T}[Fun(conj(V[:,k]).*sqrt(Σ[k]),dy) for k=1:m]

    LowRankFun(A,B)
end

## Construction in a TensorSpace via a Vector of Funs

function LowRankFun{S,T,SV}(X::Vector{Fun{S,T}},d::TensorSpace{SV,T,2})
    @assert d[1] == space(X[1])
    LowRankFun(X,d[2])
end

function LowRankFun{S,T}(X::Vector{Fun{S,T}},dy::Space)
    m=mapreduce(length,max,X)
    M=zeros(T,m,length(X))
    for k=1:length(X)
        M[1:length(X[k]),k]=X[k].coefficients
    end

    LowRankFun(M,space(X[1]),dy)
end

## Adaptive constructor selector

function LowRankFun(f::Function,dx::Space,dy::Space;method::Symbol=:standard,tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,retmax::Bool=false,gridx::Integer=64,gridy::Integer=64,maxrank::Integer=100)
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

function standardLowRankFun(f::Function,dx::Space,dy::Space;tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,gridx::Integer=64,gridy::Integer=64,maxrank::Integer=100)
    xy = checkpoints(dx⊗dy)
    T = promote_type(eltype(f(first(xy)...)),eltype(dx),eltype(domain(dx)),eltype(dy),eltype(domain(dy)))

    # We start by sampling on the given grid, find the approximate maximum and create the first rank-one approximation.
    ptsx,ptsy=points(dx,gridx),points(dy,gridy)
    X = zeros(T,gridx,gridy)
    maxabsf,r=findapproxmax!(f,X,ptsx,ptsy,gridx,gridy)
    if maxabsf < eps(zero(T))/eps(T) return LowRankFun([Fun([zero(T)],dx)],[Fun([zero(T)],dy)]),maxabsf end
    a,b=Fun(x->f(x,r[2]),dx),Fun(y->f(r[1],y),dy)

    # If necessary, we resize the grid to be at least as large as the
    # lengths of the first row and column Funs and we recompute the values of X.
    if gridx < length(a) || gridy < length(b)
        gridx,gridy = max(gridx,length(a)),max(gridy,length(b))
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

    # Eat, drink, subtract rank-one, repeat.
    for k=1:maxrank
        if norm(a.coefficients,Inf) < tol || norm(b.coefficients,Inf) < tol return LowRankFun(A,B),maxabsf end
        A,B=[A;a/sqrt(abs(a(r[1])))],[B;b/(sqrt(abs(b(r[2])))*sign(b(r[2])))]
        r=findapproxmax!(A[k],B[k],X,ptsx,ptsy,gridx,gridy)
        Ar,Br=evaluate(A,r[1]),evaluate(B,r[2])
        a,b=Fun(x->f(x,r[2]),dx,gridx) - dotu(Br,A),Fun(y->f(r[1],y),dy,gridy) - dotu(Ar,B)
        chop!(a,tol10),chop!(b,tol10)
    end
    warn("Maximum rank of " * string(maxrank) * " reached")
    return LowRankFun(A,B),maxabsf
end

## Adaptive Cholesky decomposition, when f is Hermitian positive (negative) definite

function CholeskyLowRankFun(f::Function,dx::Space;tolerance::Union{Symbol,Tuple{Symbol,Number}}=:relative,grid::Integer=64,maxrank::Integer=100)
    xy = checkpoints(dx⊗dx)
    T = promote_type(eltype(f(first(xy)...)),eltype(dx),eltype(domain(dx)))

    # We start by sampling on the given grid, find the approximate maximum and create the first rank-one approximation.
    pts=points(dx,grid)
    X = zeros(T,grid)
    maxabsf,r=findcholeskyapproxmax!(f,X,pts,grid)
    if maxabsf < eps(zero(T))/eps(T) return LowRankFun([Fun([zero(T)],dx)],[Fun([zero(T)],dy)]),maxabsf end
    a=Fun(x->f(x,r),dx)

    # If necessary, we resize the grid to be at least as large as the
    # length of the first row/column Fun and we recompute the values of X.
    if grid < length(a)
        grid = max(grid,length(a))
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

    # Eat, drink, subtract rank-one, repeat.
    for k=1:maxrank
        if norm(a.coefficients,Inf) < tol return LowRankFun(A,B),maxabsf end
        A,B=[A;a/sqrt(abs(a(r)))],[B;a/(sqrt(abs(a(r)))*sign(a(r)))]
        r=findcholeskyapproxmax!(A[k],B[k],X,pts,grid)
        Br=evaluate(B,r)
        a=Fun(x->f(x,r),dx,grid) - dotu(Br,A)
        chop!(a,tol10)
    end
    warn("Maximum rank of " * string(maxrank) * " reached")
    return LowRankFun(A,B),maxabsf
end


## Construction via TensorSpaces and ProductDomains

LowRankFun{SV,T}(f::Function,S::TensorSpace{SV,T,2};kwds...)=LowRankFun(f,S[1],S[2];kwds...)
LowRankFun(f::Function,dx::Domain,dy::Domain;kwds...)=LowRankFun(f,Space(dx),Space(dy);kwds...)
LowRankFun{D,T}(f::Function,d::ProductDomain{D,T,2};kwds...)=LowRankFun(f,d[1],d[2];kwds...)

LowRankFun(f::Function,d1::Vector,d2::Vector;kwds...)=LowRankFun(f,convert(Domain,d1),convert(Domain,d2);kwds...)
LowRankFun(f::Function;kwds...)=LowRankFun(f,Interval(),Interval();kwds...)

## Construction from values

LowRankFun{T<:Number}(A::Array{T})=LowRankFun(A,Interval{T}(),Interval{T}())
LowRankFun(c::Number,etc...)=LowRankFun((x,y)->c,etc...)

## Construction from other LowRankFuns

LowRankFun(f::LowRankFun,d1::IntervalDomain,d2::IntervalDomain)=LowRankFun(map(g->Fun(g.coefficients,d1),f.A),map(g->Fun(g.coefficients,d2),f.B))
LowRankFun(f::LowRankFun)=LowRankFun(f,Interval(),Interval())



## Utilities

function findapproxmax!(f::Function,X::Matrix,ptsx::Vector,ptsy::Vector,gridx,gridy)
    @inbounds for j=1:gridy,k=1:gridx
        X[k,j]+=f(ptsx[k],ptsy[j])
    end
    maxabsf,impt = findmax(abs(X))
    imptple = ind2sub((gridx,gridy),impt)
    maxabsf,[ptsx[imptple[1]],ptsy[imptple[2]]]
end

function findapproxmax!(A::Fun,B::Fun,X::Matrix,ptsx::Vector,ptsy::Vector,gridx,gridy)
    dX = A(ptsx)*transpose(B(ptsy))
    X[:] -= dX[:]
    maxabsf,impt = findmax(abs(X))
    imptple = ind2sub((gridx,gridy),impt)
    [ptsx[imptple[1]],ptsy[imptple[2]]]
end

function findcholeskyapproxmax!(f::Function,X::Vector,pts::Vector,grid)
    @inbounds for k=1:grid
        X[k]+=f(pts[k],pts[k])
    end
    maxabsf,impt = findmax(abs(X))
    maxabsf,pts[impt]
end

function findcholeskyapproxmax!(A::Fun,B::Fun,X::Vector,pts::Vector,grid)
    dX = A(pts).*B(pts)
    X[:] -= dX[:]
    maxabsf,impt = findmax(abs(X))
    pts[impt]
end


domain(f::LowRankFun,k::Integer)=k==1? domain(first(f.A)) : domain(first(f.B))
space(f::LowRankFun,k::Integer)=k==1? space(first(f.A)) : space(first(f.B))
space(f::LowRankFun)=f.space

Base.transpose{S,M,SS,T}(f::LowRankFun{S,M,SS,T})=LowRankFun(f.B,f.A,transpose(space(f)))

function values(f::LowRankFun)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=values(pad(f.A[k],xm))*values(pad(f.B[k],ym)).'
    end
    ret
end

#TODO: this is inconsistent with 1D where it does canonical
function coefficients(f::LowRankFun)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(f.A[k].coefficients,xm)*pad(f.B[k].coefficients,ym).'
    end
    ret
end

function coefficients(f::LowRankFun,n::Space,m::Space)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(coefficients(f.A[k],n),xm)*pad(coefficients(f.B[k],m),ym).'
    end
    ret
end

function vecpoints(f::LowRankFun,k::Integer)
    if k==1
        xm=mapreduce(length,max,f.A)
        points(space(first(f.A)),xm)
    else
        ym=mapreduce(length,max,f.B)
        points(space(first(f.B)),ym)
    end
end



evaluate{T<:Fun,M<:Fun}(A::Vector{T},B::Vector{M},x,y)=dotu(evaluate(A,x),evaluate(B,y))

evaluate(f::LowRankFun,x,y)=evaluate(f.A,f.B,x,y)
evaluate(f::LowRankFun,::Colon,::Colon)=f
evaluate(f::LowRankFun,x,::Colon)=f.B*evaluate(f.A,x)
function evaluate(f::LowRankFun,::Colon,y)
    m = maximum(map(length,f.A))
    r=rank(f)
    ret = zeros(m)

    for k=1:r
        for j=1:length(f.A[k])
            @inbounds ret[j] += f.A[k].coefficients[j]*f.B[k](y)
        end
    end

    Fun(ret,first(f.A).space)
end


## Truncate
#TODO: should reduce rank if needed
Base.chop(f::LowRankFun,tol)=LowRankFun(map(g->chop(g,tol),f.A),map(g->chop(g,tol),f.B))



## Algebra

for op = (:*,:.*,:./,:/)
    @eval ($op){T<:Fun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::LowRankFun,c::Number) = LowRankFun(($op)(f.A,c),f.B)
    @eval ($op)(c::Number,f::LowRankFun) = LowRankFun(($op)(c,f.A),f.B)
end

# Let K be a LowRankFun and f be a Fun.
# op(f,K) acts as operating in the x variable, and
# op(K,f) acts as operating in the y variable.

for op = (:*,:.*,:./,:/)
    @eval ($op){S,T,U,V}(f::Fun{S,T},A::Vector{Fun{U,V}})=map(a->($op)(f,a),A)
    @eval ($op)(f::Fun,K::LowRankFun) = LowRankFun(($op)(f,K.A),K.B)
    @eval ($op){S,T,U,V}(B::Vector{Fun{U,V}},f::Fun{S,T})=map(b->($op)(b,f),B)
    @eval ($op)(K::LowRankFun,f::Fun) = LowRankFun(K.A,($op)(K.B,f))
end

+(f::LowRankFun,g::LowRankFun)=LowRankFun([f.A,g.A],[f.B,g.B])
-(f::LowRankFun)=LowRankFun(-f.A,f.B)
-(f::LowRankFun,g::LowRankFun)=f+(-g)

## QR factorization of a LowRankFun

function Base.qr(f::LowRankFun)
    sp,r = space(f),rank(f)
    Q,R = qr(coefficients(f.A))
    BR = coefficients(f.B)*R.'
    LowRankFun(map(i->Fun(Q[:,i],sp[1]),1:r),map(i->Fun(BR[:,i],sp[2]),1:r),sp)
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
