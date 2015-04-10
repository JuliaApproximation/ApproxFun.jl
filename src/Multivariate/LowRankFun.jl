
export LowRankFun



## LowRankFun


type LowRankFun{S<:FunctionSpace,M<:FunctionSpace,T<:Number,V<:Number}<:BivariateFun
  A::Vector{Fun{S,T}}
  B::Vector{Fun{M,V}}

  function LowRankFun(A::Vector{Fun{S,T}},B::Vector{Fun{M,V}})
    @assert length(A) == length(B)
    @assert length(A) > 0
    new(A,B)
  end
end

LowRankFun{T,V,S,M}(A::Vector{Fun{S,T}},B::Vector{Fun{M,V}})=LowRankFun{S,M,T,V}(A,B)


LowRankFun{T<:Number}(A::Array{T})=LowRankFun(A,Interval(),Interval())
function LowRankFun{T<:Number,S<:FunctionSpace,W<:FunctionSpace}(X::Array{T},dx::S,dy::W)
    U,Σ,V=svd(X)
    m=max(1,count(s->s>10eps(),Σ))


    A=Fun{S,T}[Fun(U[:,k].*sqrt(Σ[k]),dx) for k=1:m]
    B=Fun{W,T}[Fun(conj(V[:,k]).*sqrt(Σ[k]),dy) for k=1:m]

    LowRankFun(A,B)
end

LowRankFun(f,S::TensorSpace)=LowRankFun(f,S[1],S[2])
LowRankFun(f,S::TensorSpace,n::Integer)=LowRankFun(f,S[1],S[2],n)
LowRankFun(f,S::TensorSpace,n::Integer,m::Integer)=LowRankFun(f,S[1],S[2],n,m)

## We take the convention that row vector pads down
# TODO: Vector pads right
for T in (:Float64,:(Complex{Float64}))
    @eval begin
        LowRankFun{S}(X::Vector{Fun{S,$T}},d::TensorSpace)=LowRankFun(X,d[1],d[2])
        function LowRankFun{S}(X::Vector{Fun{S,$T}},dy::FunctionSpace)
            m=mapreduce(length,max,X)
            M=zeros($T,m,length(X))
            for k=1:length(X)
                M[1:length(X[k]),k]=X[k].coefficients
            end

            LowRankFun(M,space(X[1]),dy)
        end
    end
end

LowRankFun(f::Function,dx::Domain,dy::Domain,nx...)=LowRankFun(f,Space(dx),Space(dy),nx...)
LowRankFun(f,d::ProductDomain)=LowRankFun(f,d[1],d[2])

function findapproxmax!(f::Function,X::Matrix,ptsx::Vector,ptsy::Vector,gridx,gridy)
    @inbounds for j=1:gridy,k=1:gridx
        X[k,j]+=f(ptsx[k],ptsy[j])
    end
    maxabsf,impt = findmax(abs(X))
    imptple = ind2sub((gridx,gridy),impt)
    maxabsf,[ptsx[imptple[1]],ptsy[imptple[2]]]
end

function LowRankFun(f::Function,dx::FunctionSpace,dy::FunctionSpace;gridx::Integer=64,gridy::Integer=64,maxrank::Integer=100)

    # We start by sampling on the given grid, find the approximate maximum and create the first rank-one approximation.
    ptsx,ptsy=points(dx,gridx),points(dy,gridy)
    X = zeros(typeof(f(ptsx[1],ptsy[1])),gridx,gridy)
    maxabsf,r=findapproxmax!(f,X,ptsx,ptsy,gridx,gridy)
    a,b=Fun(x->f(x,r[2]),dx),Fun(y->f(r[1],y),dy)

    # We resize the grid to be at least as large as the lengths of the first row and column Funs.
    gridx,gridy = max(gridx,length(a)),max(gridy,length(b))

    # If necessary, we recompute the values of X.
    if gridx == length(a) || gridy == length(b)
        ptsx,ptsy=points(dx,gridx),points(dy,gridy)
        X = zeros(typeof(f(ptsx[1],ptsy[1])),gridx,gridy)
        maxabsf,r=findapproxmax!(f,X,ptsx,ptsy,gridx,gridy)
        a,b=Fun(x->f(x,r[2]),dx),Fun(y->f(r[1],y),dy)
    end

    A,B,tol=typeof(a)[],typeof(b)[],100maxabsf*eps()

    for k=1:maxrank

        if norm(a.coefficients,Inf) < tol || norm(b.coefficients,Inf) < tol return LowRankFun(A,B) end

        A,B=[A;a/sqrt(abs(a[r[1]]))],[B;sign(b[r[2]]).*b/sqrt(abs(b[r[2]]))]

        maxabsf,r=findapproxmax!((x,y)-> - evaluate(A[end],x)*evaluate(B[end],y),X,ptsx,ptsy,gridx,gridy)

        Ar,Br=map(q->q[r[1]],A),map(q->q[r[2]],B)

        a,b=Fun(x->f(x,r[2]),dx,gridx) - dotu(Br,A),Fun(y->f(r[1],y),dy,gridy) - dotu(Ar,B)

        chop!(a,tol),chop!(b,tol)

    end
    warn("Maximum rank of " * string(maxrank) * " reached")
    return LowRankFun(A,B)
end



LowRankFun(f::Function,d1::Vector,d2::Vector)=LowRankFun(f,Interval(d1),Interval(d2))
LowRankFun(f::Function)=LowRankFun(f,Interval(),Interval())

LowRankFun(f::LowRankFun,d1::IntervalDomain,d2::IntervalDomain)=LowRankFun(map(g->Fun(g.coefficients,d1),f.A),map(g->Fun(g.coefficients,d2),f.B))

LowRankFun(f::LowRankFun)=LowRankFun(f,Interval(),Interval())


domain(f::LowRankFun,k::Integer)=k==1? domain(first(f.A)) : domain(first(f.B))
space(f::LowRankFun,k::Integer)=k==1? space(first(f.A)) : space(first(f.B))


Base.size(f::LowRankFun,k::Integer)=k==1?mapreduce(length,max,f.A):mapreduce(length,max,f.B)
Base.size(f::LowRankFun)=size(f,1),size(f,2)

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

function coefficients(f::LowRankFun,n::FunctionSpace,m::FunctionSpace)
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




evaluate(f::LowRankFun,x::Real,y::Real)=evaluate(f.A,f.B,x,y)
evaluate(f::LowRankFun,x::Real,::Colon)=f.B*evaluate(f.A,x)
function evaluate(f::LowRankFun,::Colon,y::Real)
    m = maximum(map(length,f.A))
    r=rank(f)
    ret = zeros(m)

    for k=1:r
        for j=1:length(f.A[k])
            @inbounds ret[j] += f.A[k].coefficients[j]*f.B[k][y]
        end
    end

    Fun(ret,first(f.A).space)
end



Base.rank(f::LowRankFun)=length(f.A)
evaluate{T<:Fun,M<:Fun}(A::Vector{T},B::Vector{M},x,y)=dotu(evaluate(A,x),evaluate(B,y))

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


## Special functions

Base.real(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(real,u.B),-map(imag,u.B)])
Base.imag(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(imag,u.B),map(real,u.B)])


## Calculus


Base.sum(g::LowRankFun)=dotu(map(sum,g.A),map(sum,g.B))
Base.sum(g::LowRankFun,n::Integer)=(n==1)?dotu(g.B,map(sum,g.A)):dotu(g.A,map(sum,g.B))
Base.cumsum(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(cumsum,g.A),copy(g.B)):LowRankFun(copy(g.A),map(cumsum,g.B))
integrate(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(integrate,g.A),copy(g.B)):LowRankFun(copy(g.A),map(integrate,g.B))


