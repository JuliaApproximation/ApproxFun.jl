
# represents function of the form r^m g(r^2)
# as r^m P^{.5+a,b-.5}(r^2)

immutable JacobiSquareSpace <: IntervalDomainSpace
    m::Int
    a::Int
    b::Int
    domain::Union(Interval,AnyDomain)
end
JacobiSquareSpace(m::Integer,d::Domain)=JacobiSquareSpace(m,m,0,d)
JacobiSquareSpace(m::Integer)=JacobiSquareSpace(m,Interval(1.,0.))


jacobispace(B::JacobiSquareSpace)=JacobiSpace(B.a+.5,B.b-.5,domain(B))


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


function addentries!{T}(M::Multiplication{JacobiWeightSpace{ChebyshevSpace},JacobiSquareSpace,T},A::ShiftArray,kr::Range)
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    addentries!(ConstantOperator(0.5M.f.coefficients[1]),A,kr)
end
function rangespace{T}(M::Multiplication{JacobiWeightSpace{ChebyshevSpace},JacobiSquareSpace,T})
    @assert length(M.f)==1
    @assert M.f.space.α ==0.
    @assert isinteger(M.f.space.β)
    ds=domainspace(M)
    JacobiSquareSpace(ds.m+int(M.f.space.β),ds.a,ds.b,domain(M))
end




function Derivative(S::JacobiSquareSpace)
     # we have D[r^m f(r^2)] = r^{m-1} (m f(r^2) + 2r^2 f'(r^2))
     a=S.a;b=S.b;m=S.m
     d=domain(S)
     @assert d==Interval(1.,0.)
     
     JS=jacobispace(S)
     D=Derivative(JS)
     M=Multiplication(Fun(identity,d),rangespace(D))
     
     DerivativeWrapper(SpaceOperator(m+2*M*D,S,JacobiSquareSpace(m-1,a+1,b+1,d)),1)
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


##TODO:ConversionWrapper
function addentries!{T}(C::Conversion{JacobiSquareSpace,JacobiSquareSpace,T},SA::ShiftArray,kr::Range)
    dm=domain(C)
    A=domainspace(C);B=rangespace(C)
    
    addentries!(Conversion(jacobispace(A),jacobispace(B)),SA,kr)
end

function bandinds{T}(C::Conversion{JacobiSquareSpace,JacobiSquareSpace,T})
    dm=domain(C)
    A=domainspace(C);B=rangespace(C)
    
    bandinds(Conversion(jacobispace(A),jacobispace(B)))
end



function Base.getindex(op::Evaluation{JacobiSquareSpace,Bool},kr::Range)
    @assert !op.x && op.order ==0
    getindex(Evaluation(jacobispace(op.space),false,0),kr)
end
