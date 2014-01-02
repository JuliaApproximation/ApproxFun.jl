
export coefficientmatrix, Fun2D


## Vector of fun routines

function coefficientmatrix{T<:IFun}(B::Vector{T})
    m=mapreduce(length,max,B)
  n=length(B)
  ret = Array(Float64,m,length(B))
  for k=1:m,j=1:n
    ret[k,j] = length(B[j]) < k ? 0. : B[j].coefficients[k]
  end
  
  ret
end

*{T<:IFun}(v::Vector{T},a::Vector)=IFun(coefficientmatrix(v)*a,first(v).domain)



## Fun2D


type Fun2D{T<:IFun}
#TODO: allow mix of IFun and FFun
  A::Vector{T}
  B::Vector{T}
  
  function Fun2D(A::Vector{T},B::Vector{T})
    @assert length(A) == length(B)
    new(A,B)
  end
end

Fun2D{T<:IFun}(A::Vector{T},B::Vector{T})=Fun2D{T}(A,B)



  
function findapproxmax(f::Function,dx::IntervalDomain,dy::IntervalDomain)
ptsx=points(dx,20)
ptsy=points(dy,20)

  mpt=[fromcanonical(dx,rand()),fromcanonical(dy,rand())]  
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

function Fun2D(f::Function,dx::IntervalDomain,dy::IntervalDomain)
    tol=1000eps()
    
    r=findapproxmax(f,dx,dy)
    a=Fun(x->f(x,r[2]),dx)
    b=Fun(y->f(r[1],y),dy)
    A=typeof(a)[];B=typeof(b)[];
    
    
    while norm(a) > tol || norm(b) > tol
        A=[A,a];B=[B,b/b[r[2]]]    
        r=findapproxmax((x,y)->f(x,y) - evaluate(A,B,x,y),dx,dy)
        Ar=map(q->q[r[1]],A)
        Br=map(q->q[r[2]],B)
        a=Fun(x->f(x,r[2]),dx) - A*Br
        b=Fun(y->f(r[1],y),dy)- B*Ar
    end
      
    Fun2D(A,B)
end

Fun2D(f::Function,d1::Vector,d2::Vector)=Fun2D(f,Interval(d1),Interval(d2))
Fun2D(f::Function)=Fun2D(f,Interval(),Interval())

domain(f::Fun2D,k::Integer)=k==1? first(f.A).domain : first(f.B).domain


evaluate{T<:IFun}(A::Vector{T},B::Vector{T},x::Real,y::Real)=dot(evaluate(A,x),evaluate(B,y))
evaluate{T<:IFun}(A::Vector{T},x::Real)=Float64[A[k][x] for k=1:length(A)]
function evaluate{T<:IFun}(A::Vector{T},x::Vector{Float64})
    x = tocanonical(first(A),x)

    n=length(x)
    ret=Array(Float64,length(A),n)
    
    bk=Array(Float64,n)
    bk1=Array(Float64,n)
    bk2=Array(Float64,n)
    
    for k=1:length(A)
        bk=clenshaw(A[k].coefficients,x,bk,bk1,bk2)
        
        for j=1:n
            ret[k,j]=bk[j]
        end
    end
    
    ret
end

evaluate(f::Fun2D,x::Real,y::Real)=evaluate(f.A,f.B,x,y)
evaluate(f::Fun2D,x::Real,::Colon)=dot(evaluate(f.A,x),f.B)
function evaluate(f::Fun2D,::Colon,y::Real)
    m = maximum(map(length,f.A))
    r=rank(f)
    ret = zeros(m)
    ret_v = unsafe_view(ret)
    
    for k=1:r
        for j=1:length(f.A[k])
            ret_v[j] += f.A[k].coefficients[j]*f.B[k][y]
        end
    end
    
    IFun(ret,first(f.A).domain)
end

Base.getindex(f::Fun2D,x,y)=evaluate(f,x,y)

Base.rank(f::Fun2D)=length(f.A)
Base.sum(g::Fun2D)=dot(map(sum,g.A),map(sum,g.B))
Base.sum(g::Fun2D,n::Integer)=(n==1)?dot(map(sum,g.A),g.B):dot(map(sum,g.B),g.A)
Base.cumsum(g::Fun2D,n::Integer)=(n==1)?Fun2D(map(cumsum,g.A),copy(g.B)):Fun2D(copy(g.A),map(cumsum,g.B))

for op = (:*,:.*,:./,:/)
    @eval ($op){T<:IFun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::Fun2D,c::Number) = Fun2D(($op)(f.A,c),f.B)
end 


