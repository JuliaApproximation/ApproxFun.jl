
# represents function of the form r^m g(r^2)
# as r^m P^{.5+a,b-.5}(r^2)

immutable JacobiSquare <: IntervalSpace
    m::Int
    a::Int
    b::Int
    domain::Union(Interval,AnyDomain)
end
JacobiSquare(m::Integer,d::Domain)=JacobiSquare(m,m,0,d)
JacobiSquare(m::Integer)=JacobiSquare(m,Interval(1.,0.))


jacobispace(B::JacobiSquare)=Jacobi(B.a+.5,B.b-.5,domain(B))


spacescompatible(A::JacobiSquare,B::JacobiSquare)=A.a==B.a&&A.m==B.m&&A.b==B.b

canonicalspace(S::JacobiSquare)=S


# We assume domains is [1.,0.]

points(S::JacobiSquare,n)=sqrt(fromcanonical(S.domain,gausschebyshev(n,4)[1]))

plan_transform(S::JacobiSquare,n)=gausschebyshev(n,4)

transform(S::JacobiSquare,vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))


function transform(S::JacobiSquare,vals::Vector,xw::@compat(Tuple{Vector,Vector}))
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
        V=jacobip(0:n-div(m,2)-1,m+0.5,-0.5,x)'
        nrm=(V.^2)*(w2.*mw)
        (V*(mw.*vals))./nrm*2^(m/2)
    else
        error("transform only implemented for first case")
    end
end


#evaluate{m,order}(f::Fun{JacobiSquare{m,order}},x::Number)=x^m*dot(jacobip(0:length(f)-1,m+0.5+order,order-0.5,tocanonical(f,x^2)),f.coefficients)*2^(m/2)



plan_itransform(S::JacobiSquare,n)=points(S,n)
#TODO: general domain
itransform(S::JacobiSquare,cfs::Vector)=itransform(S,cfs,plan_itransform(S,length(cfs)))
itransform(S::JacobiSquare,cfs::Vector,x)=isempty(cfs)?zeros(x):x.^S.m.*jacobip(0:length(cfs)-1,S.a+0.5,S.b-0.5,tocanonical(S,x.^2))*cfs
evaluate{T}(f::Fun{JacobiSquare,T},x::Vector)=itransform(f.space,f.coefficients,x)
evaluate{T}(f::Fun{JacobiSquare,T},x::Number)=itransform(f.space,f.coefficients,[x])[1]



## Operators

# Override JacobiWeight default
Multiplication{T}(f::Fun{JacobiWeight{Chebyshev},T},S::JacobiSquare)=Multiplication{JacobiWeight{Chebyshev},JacobiSquare,T,T}(f,S)
bandinds{T}(M::Multiplication{JacobiWeight{Chebyshev},JacobiSquare,T})=0,0

function addentries!{T}(M::Multiplication{JacobiWeight{Chebyshev},JacobiSquare,T},A,kr::Range)
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    addentries!(ConstantOperator(2.0^M.f.space.β*M.f.coefficients[1]),A,kr)
end
function rangespace{T}(M::Multiplication{JacobiWeight{Chebyshev},JacobiSquare,T})
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    @assert isinteger(M.f.space.β)
    ds=domainspace(M)
    JacobiSquare(ds.m+round(Int,M.f.space.β),ds.a,ds.b,domain(M))
end




function Derivative(S::JacobiSquare)

     a=S.a;b=S.b;m=S.m
     d=domain(S)
     @assert d==Interval(1.,0.)

     JS=jacobispace(S)
     D=Derivative(JS)
     M=Multiplication(2Fun(identity,d),rangespace(D))

      if m==0
      # we have D[ f(r^2)] = 2r f'(r^2) = 2 r^(-1)*r^2 f'(r^2)
         DerivativeWrapper(SpaceOperator(M*D,S,JacobiSquare(-1,a+1,b+1,d)),1)
      else
     # we have D[r^m f(r^2)] = r^{m-1} (m f(r^2) + 2r^2 f'(r^2))
        DerivativeWrapper(SpaceOperator(m+M*D,S,JacobiSquare(m-1,a+1,b+1,d)),1)
    end
end


function Derivative(S::JacobiSquare,k::Integer)
     if k==1
         Derivative(S)
     else
         D=Derivative(S)
         DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
     end
end



# return the space that has banded Conversion to the other
conversion_rule(A::JacobiSquare,B::JacobiSquare)=JacobiSquare(max(A.m,B.m),min(A.a,B.a),min(A.b,B.b),domain(A))
maxspace(A::JacobiSquare,B::JacobiSquare)=JacobiSquare(min(A.m,B.m),max(A.a,B.a),max(A.b,B.b),domain(A))


##TODO:ConversionWrapper

function Conversion(A::JacobiSquare,B::JacobiSquare)
    if A.m==B.m
        ConversionWrapper(SpaceOperator(Conversion(jacobispace(A),jacobispace(B)),A,B))
    else
        @assert A.m > B.m && iseven(A.m-B.m)
        r=Fun(identity,domain(B))
        M=Multiplication(r.^div(A.m-B.m,2),jacobispace(B)) #this is multiplication by r^(2*p)
        ConversionWrapper(SpaceOperator(M*Conversion(jacobispace(A),jacobispace(B)),A,B))
    end
end




function Base.getindex(op::Evaluation{JacobiSquare,Bool},kr::Range)
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
