
export LowRankFun



## LowRankFun


type LowRankFun{T<:Number,S<:DomainSpace,M<:DomainSpace}<:MultivariateFun
  A::Vector{Fun{S,T}}
  B::Vector{Fun{M,T}}
  
  function LowRankFun(A::Vector{Fun{S,T}},B::Vector{Fun{M,T}})
    @assert length(A) == length(B)
    @assert length(A) > 0
    new(A,B)
  end
end

LowRankFun{T,S,M}(A::Vector{Fun{S,T}},B::Vector{Fun{M,T}})=LowRankFun{T,S,M}(A,B)


LowRankFun{T<:Number}(A::Array{T})=LowRankFun(A,Interval(),Interval())
function LowRankFun{T<:Number,S<:FunctionSpace,W<:FunctionSpace}(X::Array{T},dx::S,dy::W)
    U,Σ,V=svd(X)
    m=max(1,count(s->s>10eps(),Σ))
    

    A=Fun{S,T}[Fun(U[:,k].*sqrt(Σ[k]),dx) for k=1:m]
    B=Fun{W,T}[Fun(conj(V[:,k]).*sqrt(Σ[k]),dy) for k=1:m]

    LowRankFun(A,B)
end

## We take the convention that row vector pads down
# TODO: Vector pads right
for T in (:Float64,:(Complex{Float64}))
    @eval begin
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


findapproxmax(f::Function,dx::FunctionSpace,dy::FunctionSpace)=findapproxmax(f,dx,dy, 40, 40)
function findapproxmax(f::Function,dx::FunctionSpace,dy::FunctionSpace, gridx::Integer, gridy::Integer)
    ptsx=points(dx,gridx)
    ptsy=points(dy,gridy)

  mpt=[fromcanonical(dx,0.),fromcanonical(dy,0.)]  
  maxi=abs(f(mpt[1],mpt[2]))


    for k = 1:length(ptsx),j=1:length(ptsy)
      val=abs(f(ptsx[k],ptsy[j])) 
      if val > maxi
        maxi = val
        mpt[1]=ptsx[k];mpt[2]=ptsy[j]
      end
    end
    mpt
end
LowRankFun(f::Function,dx::FunctionSpace,dy::FunctionSpace)=LowRankFun(f,dx,dy,40,40)

LowRankFun(f::Function,dx::Domain,dy::Domain,nx...)=LowRankFun(f,Space(dx),Space(dy),nx...)
LowRankFun(f,d::ProductDomain)=LowRankFun(f,d[1],d[2])
function LowRankFun(f::Function,dx::FunctionSpace,dy::FunctionSpace,gridx::Integer,gridy::Integer;maxrank=100::Integer)
    tol=1000eps()
    
    r=findapproxmax(f,dx,dy,gridx,gridy)
    a=Fun(x->f(x,r[2]),dx)
    b=Fun(y->f(r[1],y),dy)
    A=typeof(a)[];B=typeof(b)[];
    
    
    for k=1:maxrank
        if norm(a.coefficients) < tol && norm(b.coefficients) < tol
            return LowRankFun(A,B)
        end
        
        
        ##Todo negative orientation 
        A=[A,a/sqrt(abs(a[r[1]]))];B=[B,sign(b[r[2]]).*b/sqrt(abs(b[r[2]]))]    
        r=findapproxmax((x,y)->f(x,y) - evaluate(A,B,x,y),dx,dy,gridx,gridy)
        Ar=map(q->q[r[1]],A)
        Br=map(q->q[r[2]],B)
        
        ##TODO: Should allow FFun
        a=Fun(x->f(x,r[2]),dx; method="abszerocoefficients") - dot(conj(Br),A)
        b=Fun(y->f(r[1],y),dy; method="abszerocoefficients")- dot(conj(Ar),B)
        
        
        ##Remove coefficients that get killed by a/b
        a=chop!(a,10*sqrt(abs(a[r[1]]))*eps()/maximum(abs(b.coefficients)))
        b=chop!(b,10*sqrt(abs(b[r[2]]))*eps()/maximum(abs(a.coefficients)))        
    end
      
    error("Maximum rank of " * string(maxrank) * " reached")
end

LowRankFun(f::Function,d1::Vector,d2::Vector)=LowRankFun(f,Interval(d1),Interval(d2))
LowRankFun(f::Function)=LowRankFun(f,Interval(),Interval())

LowRankFun(f::LowRankFun,d1::IntervalDomain,d2::IntervalDomain)=LowRankFun(map(g->Fun(g.coefficients,d1),f.A),map(g->Fun(g.coefficients,d2),f.B))

LowRankFun(f::LowRankFun)=LowRankFun(f,Interval(),Interval())


domain(f::LowRankFun,k::Integer)=k==1? domain(first(f.A)) : domain(first(f.B))
space(f::LowRankFun,k::Integer)=k==1? space(first(f.A)) : space(first(f.B))


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

function coefficients(f::LowRankFun,n::Integer,m::Integer)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)    
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(coefficients(f.A[k],n),xm)*pad(coefficients(f.B[k],m),ym).'
    end
    ret
end

function points(f::LowRankFun,k::Integer)
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
Base.sum(g::LowRankFun)=dot(map(sum,g.A),map(sum,g.B)) #TODO: not complexconjugate
evaluate{T<:Fun,M<:Fun}(A::Vector{T},B::Vector{M},x,y)=dot(evaluate(A,x),evaluate(B,y)) #TODO: not complexconjugate


Base.sum(g::LowRankFun,n::Integer)=(n==1)?dotu(g.B,map(sum,g.A)):dotu(g.A,map(sum,g.B))
Base.cumsum(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(cumsum,g.A),copy(g.B)):LowRankFun(copy(g.A),map(cumsum,g.B))
integrate(g::LowRankFun,n::Integer)=(n==1)?LowRankFun(map(integrate,g.A),copy(g.B)):LowRankFun(copy(g.A),map(integrate,g.B))

for op = (:*,:.*,:./,:/)
    @eval ($op){T<:Fun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::LowRankFun,c::Number) = LowRankFun(($op)(f.A,c),f.B)
    @eval ($op)(c::Number,f::LowRankFun) = LowRankFun(($op)(c,f.A),f.B)
end 

+(f::LowRankFun,g::LowRankFun)=LowRankFun([f.A,g.A],[f.B,g.B])
-(f::LowRankFun)=LowRankFun(-f.A,f.B)
-(f::LowRankFun,g::LowRankFun)=f+(-g)


Base.real(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(real,u.B),-map(imag,u.B)])
Base.imag(u::LowRankFun)=LowRankFun([map(real,u.A),map(imag,u.A)],[map(imag,u.B),map(real,u.B)])
