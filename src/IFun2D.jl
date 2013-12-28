


type IFun2D
  A
  B
end



  

function IFun2D(f::Function)
    tol=1000eps()
    
    r=rand(2)
    a=Fun(x->f(x,r[2]))
    b=Fun(y->f(r[1],y))
    A=[];B=[];
    
    
    while norm(a) > tol || norm(b) > tol
        A=[A,a];B=[B,b/b[r[2]]]    
        r=rand(2)
        Ar=map(q->q[r[1]],A)
        Br=map(q->q[r[2]],B)
        a=Fun(x->f(x,r[2])) - dot(Br,A)
        b=Fun(y->f(r[1],y))-dot(Ar,B)
    end
      
    IFun2D(A,B)
end

Base.getindex(f::IFun2D,x,y)=dot(map(q->q[x],f.A),map(q->q[y],f.B))


