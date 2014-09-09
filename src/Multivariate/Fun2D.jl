
export Fun2D



## Fun2D


type Fun2D{T<:AbstractFun,M<:AbstractFun}<:MultivariateFun
  A::Vector{T}
  B::Vector{M}
  
  function Fun2D(A::Vector{T},B::Vector{M})
    @assert length(A) == length(B)
    @assert length(A) > 0
    new(A,B)
  end
end

Fun2D{T<:AbstractFun,M<:AbstractFun}(A::Vector{T},B::Vector{M})=Fun2D{T,M}(A,B)


Fun2D{T<:Number}(A::Array{T})=Fun2D(A,Interval(),Interval())
function Fun2D{T<:Number}(X::Array{T},dx::IntervalDomain,dy::IntervalDomain)
    U,Σ,V=svd(X)
    m=max(1,count(s->s>10eps(),Σ))
    

    A=IFun[IFun(U[:,k].*sqrt(Σ[k]),dx) for k=1:m]
    B=IFun[IFun(conj(V[:,k]).*sqrt(Σ[k]),dy) for k=1:m]

    Fun2D(A,B)
end

## We take the convention that row vector pads down
# TODO: Vector pads right
function Fun2D{T<:Number}(X::Array{IFun{T},2},dy::IntervalDomain)
    @assert size(X,1)==1
    
    m=mapreduce(length,max,X)
    M=zeros(T,m,length(X))
    for k=1:length(X)
        M[1:length(X[k]),k]=X[k].coefficients
    end
    
    Fun2D(M,domain(X[1]),dy)
end


findapproxmax(f::Function,dx::Domain,dy::Domain)=findapproxmax(f,dx,dy, 40, 40)
function findapproxmax(f::Function,dx::Domain,dy::Domain, gridx::Integer, gridy::Integer)
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


Fun2D(f::Function,dx::Domain,dy::Domain)=Fun2D(f,dx,dy,40,40)
function Fun2D(f::Function,dx::Domain,dy::Domain,gridx::Integer,gridy::Integer;maxrank=100::Integer)
    tol=1000eps()
    
    r=findapproxmax(f,dx,dy,gridx,gridy)
    a=Fun(x->f(x,r[2]),dx)
    b=Fun(y->f(r[1],y),dy)
    A=typeof(a)[];B=typeof(b)[];
    
    
    for k=1:maxrank
        if norm(a.coefficients) < tol && norm(b.coefficients) < tol
            return Fun2D(A,B)
        end
        
        
        ##Todo negative orientation 
        A=[A,a/sqrt(abs(a[r[1]]))];B=[B,sign(b[r[2]]).*b/sqrt(abs(b[r[2]]))]    
        r=findapproxmax((x,y)->f(x,y) - evaluate(A,B,x,y),dx,dy,gridx,gridy)
        Ar=map(q->q[r[1]],A)
        Br=map(q->q[r[2]],B)
        
        ##TODO: Should allow FFun
        a=IFun(x->f(x,r[2]),dx; method="abszerocoefficients") - A*Br
        b=IFun(y->f(r[1],y),dy; method="abszerocoefficients")- B*Ar
        
        
        ##Remove coefficients that get killed by a/b
        a=chop!(a,10*sqrt(abs(a[r[1]]))*eps()/maximum(abs(b.coefficients)))
        b=chop!(b,10*sqrt(abs(b[r[2]]))*eps()/maximum(abs(a.coefficients)))        
    end
      
    error("Maximum rank of " * string(maxrank) * " reached")
end

Fun2D(f::Function,d1::Vector,d2::Vector)=Fun2D(f,Interval(d1),Interval(d2))
Fun2D(f::Function)=Fun2D(f,Interval(),Interval())

Fun2D(f::Fun2D,d1::IntervalDomain,d2::IntervalDomain)=Fun2D(map(g->Fun(g.coefficients,d1),f.A),map(g->Fun(g.coefficients,d2),f.B))

Fun2D(f::Fun2D)=Fun2D(f,Interval(),Interval())


domain(f::Fun2D,k::Integer)=k==1? domain(first(f.A)) : domain(first(f.B))

function values(f::Fun2D)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)    
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=values(pad(f.A[k],xm))*values(pad(f.B[k],ym)).'
    end
    ret
end

function coefficients(f::Fun2D)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)    
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(f.A[k].coefficients,xm)*pad(f.B[k].coefficients,ym).'
    end
    ret
end

function coefficients(f::Fun2D,n::Integer,m::Integer)
    xm=mapreduce(length,max,f.A)
    ym=mapreduce(length,max,f.B)    
    ret=zeros(xm,ym)
    for k=1:length(f.A)
        ret+=pad(coefficients(f.A[k],n),xm)*pad(coefficients(f.B[k],m),ym).'
    end
    ret
end

function points(f::Fun2D,k::Integer)
    if k==1
        xm=mapreduce(length,max,f.A)
        points(first(f.A),xm)
    else
        ym=mapreduce(length,max,f.B)
        points(first(f.B),ym)
    end
end



evaluate(f::Fun2D,x::Real,y::Real)=evaluate(f.A,f.B,x,y)
evaluate(f::Fun2D,x::Real,::Colon)=f.B*evaluate(f.A,x)
function evaluate(f::Fun2D,::Colon,y::Real)
    m = maximum(map(length,f.A))
    r=rank(f)
    ret = zeros(m)
    
    for k=1:r
        for j=1:length(f.A[k])
            @inbounds ret[j] += f.A[k].coefficients[j]*f.B[k][y]
        end
    end
    
    IFun(ret,first(f.A).space)
end



Base.rank(f::Fun2D)=length(f.A)
Base.sum(g::Fun2D)=dot(map(sum,g.A),map(sum,g.B)) #TODO: not complexconjugate
evaluate{T<:AbstractFun,M<:AbstractFun}(A::Vector{T},B::Vector{M},x,y)=dot(evaluate(A,x),evaluate(B,y)) #TODO: not complexconjugate


Base.sum(g::Fun2D,n::Integer)=(n==1)?g.B*map(sum,g.A):g.A*map(sum,g.B) #TODO: Fun*vec should be Array[IFun]
Base.cumsum(g::Fun2D,n::Integer)=(n==1)?Fun2D(map(cumsum,g.A),copy(g.B)):Fun2D(copy(g.A),map(cumsum,g.B))
integrate(g::Fun2D,n::Integer)=(n==1)?Fun2D(map(integrate,g.A),copy(g.B)):Fun2D(copy(g.A),map(integrate,g.B))

for op = (:*,:.*,:./,:/)
    @eval ($op){T<:IFun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::Fun2D,c::Number) = Fun2D(($op)(f.A,c),f.B)
    @eval ($op)(c::Number,f::Fun2D) = Fun2D(($op)(c,f.A),f.B)
end 

+(f::Fun2D,g::Fun2D)=Fun2D([f.A,g.A],[f.B,g.B])
-(f::Fun2D)=Fun2D(-f.A,f.B)
-(f::Fun2D,g::Fun2D)=f+(-g)


real(u::Fun2D)=Fun2D([map(real,u.A),map(imag,u.A)],[map(real,u.B),-map(imag,u.B)])
imag(u::Fun2D)=Fun2D([map(real,u.A),map(imag,u.A)],[map(imag,u.B),map(real,u.B)])
