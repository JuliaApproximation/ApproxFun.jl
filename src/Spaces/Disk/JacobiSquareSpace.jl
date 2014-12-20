
# represents function of the form r^m g(r^2)
# as r^m P^{.5+a,b-.5}(r^2)

immutable JacobiSquareSpace <: IntervalSpace
    m::Int
    a::Int
    b::Int
    domain::Union(Interval,AnyDomain)
end
JacobiSquareSpace(m::Integer,d::Domain)=JacobiSquareSpace(m,m,0,d)
JacobiSquareSpace(m::Integer)=JacobiSquareSpace(m,Interval(1.,0.))


jacobispace(B::JacobiSquareSpace)=Jacobi(B.a+.5,B.b-.5,domain(B))


spacescompatible(A::JacobiSquareSpace,B::JacobiSquareSpace)=A.a==B.a&&A.m==B.m&&A.b==B.b

# We assume domains is [1.,0.]

points(S::JacobiSquareSpace,n)=sqrt(fromcanonical(S.domain,gausschebyshev(n,4)[1]))

plan_transform(S::JacobiSquareSpace,n)=gausschebyshev(n,4)

transform(S::JacobiSquareSpace,vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))


function transform(S::JacobiSquareSpace,vals::Vector,xw::(Vector,Vector))
    m=S.m;a=S.a;b=S.b
    x,w=xw
        
    n=length(vals)
    
        
    if m==a==b==0

        V=jacobip(0:n-1,0.5,-0.5,x)'
        nrm=(V.^2)*w
        (V*(w.*vals))./nrm
    elseif m==a && b==0
        ## Same as jacobitransform.jl
    
        w2=(1-x).^(m/2)
        mw=w2.*w
        V=jacobip(0:n-int(m/2)-1,m+0.5,-0.5,x)'  
        nrm=(V.^2)*(w2.*mw)    
        (V*(mw.*vals))./nrm*2^(m/2)    
    else
        error("transform only implemented for first case")
    end
end


#evaluate{m,order}(f::Fun{JacobiSquareSpace{m,order}},x::Number)=x^m*dot(jacobip(0:length(f)-1,m+0.5+order,order-0.5,tocanonical(f,x^2)),f.coefficients)*2^(m/2)



plan_itransform(S::JacobiSquareSpace,n)=points(S,n)
#TODO: general domain
itransform(S::JacobiSquareSpace,cfs::Vector)=itransform(S,cfs,plan_itransform(S,length(cfs)))
itransform(S::JacobiSquareSpace,cfs::Vector,x)=x.^S.m.*jacobip(0:length(cfs)-1,S.a+0.5,S.b-0.5,tocanonical(S,x.^2))*cfs
evaluate{T}(f::Fun{JacobiSquareSpace,T},x::Vector)=itransform(f.space,f.coefficients,x)
evaluate{T}(f::Fun{JacobiSquareSpace,T},x::Number)=itransform(f.space,f.coefficients,[x])[1]



## Operators


function addentries!{T}(M::Multiplication{JacobiWeight{Chebyshev},JacobiSquareSpace,T},A::ShiftArray,kr::Range)
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    addentries!(ConstantOperator(2.0^M.f.space.β*M.f.coefficients[1]),A,kr)
end
function rangespace{T}(M::Multiplication{JacobiWeight{Chebyshev},JacobiSquareSpace,T})
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    @assert isinteger(M.f.space.β)
    ds=domainspace(M)
    JacobiSquareSpace(ds.m+int(M.f.space.β),ds.a,ds.b,domain(M))
end




function Derivative(S::JacobiSquareSpace)

     a=S.a;b=S.b;m=S.m
     d=domain(S)
     @assert d==Interval(1.,0.)
     
     JS=jacobispace(S)
     D=Derivative(JS)
     
      if m==0
      # we have D[ f(r^2)] = 2r f'(r^2)
         DerivativeWrapper(SpaceOperator(2*D,S,JacobiSquareSpace(1,a+1,b+1,d)),1)     
      else
     # we have D[r^m f(r^2)] = r^{m-1} (m f(r^2) + 2r^2 f'(r^2))     
        M=Multiplication(Fun(identity,d),rangespace(D))
        DerivativeWrapper(SpaceOperator(m+2*M*D,S,JacobiSquareSpace(m-1,a+1,b+1,d)),1)
    end
end
     

function Derivative(S::JacobiSquareSpace,k::Integer)
     if k==1
         Derivative(S)
     else
         D=Derivative(S)
         DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
     end
end



# return the space that has banded Conversion to the other
function conversion_rule(A::JacobiSquareSpace,B::JacobiSquareSpace)
    if A.m>=B.m && A.a<= B.a && A.b<= B.b
        A
    else
        NoSpace()
    end
end

maxspace(A::JacobiSquareSpace,B::JacobiSquareSpace)=JacobiSquareSpace(min(A.m,B.m),max(A.a,B.a),max(A.b,B.b),domain(A))
minspace(A::JacobiSquareSpace,B::JacobiSquareSpace)=JacobiSquareSpace(max(A.m,B.m),min(A.a,B.a),min(A.b,B.b),domain(A))


##TODO:ConversionWrapper

function Conversion(A::JacobiSquareSpace,B::JacobiSquareSpace)
    if A.m==B.m
        ConversionWrapper(SpaceOperator(Conversion(jacobispace(A),jacobispace(B)),A,B))
    else
        @assert A.m > B.m && iseven(A.m-B.m)
        r=Fun(identity,domain(B))
        M=Multiplication(r.^div(A.m-B.m,2),jacobispace(B)) #this is multiplication by r^(2*p)
        ConversionWrapper(SpaceOperator(M*Conversion(jacobispace(A),jacobispace(B)),A,B))        
    end
end




function Base.getindex(op::Evaluation{JacobiSquareSpace,Bool},kr::Range)
    @assert !op.x && op.order <= 1
    m=op.space.m
    js=jacobispace(op.space)
    if op.order ==0
        getindex(Evaluation(js,false,0),kr)
    elseif m==0
        2getindex(Evaluation(js,false,1),kr)
    else
        2getindex(Evaluation(js,false,1),kr)+m*getindex(Evaluation(js,false,0),kr)
    end
end
