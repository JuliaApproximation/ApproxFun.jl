
export Fun2D

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



  

function Fun2D(f::Function)
    tol=1000eps()
    
    r=rand(2)
    a=Fun(x->f(x,r[2]))
    b=Fun(y->f(r[1],y))
    A=typeof(a)[];B=typeof(b)[];
    
    
    while norm(a) > tol || norm(b) > tol
        A=[A,a];B=[B,b/b[r[2]]]    
        r=rand(2)
        Ar=map(q->q[r[1]],A)
        Br=map(q->q[r[2]],B)
        a=Fun(x->f(x,r[2])) - dot(Br,A)
        b=Fun(y->f(r[1],y))-dot(Ar,B)
    end
      
    Fun2D(A,B)
end

evaluate(f::Fun2D,x::Real,y::Real)=dot(map(q->q[x],f.A),map(q->q[y],f.B))
evaluate(f::Fun2D,x::Real,::Colon)=dot(map(q->q[x],f.A),f.B)
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
    
    IFun(ret,f.A[1].domain)
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

function sample(f::Fun2D,n::Integer)
    ry=sample(sum(f,1),n)
    rx=[sample(evaluate(f,:,ry[k])) for k=1:n]
    [rx ry]
end

sample(f::Fun2D)=sample(f,1)[1,:]
