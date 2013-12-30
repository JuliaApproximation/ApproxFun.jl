
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

Base.getindex(f::Fun2D,x,y)=dot(map(q->q[x],f.A),map(q->q[y],f.B))
Base.rank(f::Fun2D)=length(f.A)
Base.sum(g::Fun2D)=dot(map(sum,g.A),map(sum,g.B))

for op = (:*,:.*,:./,:/)
    @eval ($op){T<:IFun}(A::Array{T,1},c::Number)=map(f->($op)(f,c),A)
    @eval ($op)(f::Fun2D,c::Number) = Fun2D(($op)(f.A,c),f.B)
end 
