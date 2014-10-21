
# represents function of the form r^m g(r^2)
# as r^m P^{.5+a,b-.5}(r^2)

immutable JacobiSquareSpace{m,a,b} <: IntervalDomainSpace
    domain::Union(Interval,AnyDomain)
end
JacobiSquareSpace(m,d::Domain)=JacobiSquareSpace{m,m,0}(d)
JacobiSquareSpace(m::Integer)=JacobiSquareSpace(m,Interval(1.,0.))
spacescompatible{m,a,b}(::JacobiSquareSpace{m,a,b},::JacobiSquareSpace{m,a,b})=true

# We assume domains is [1.,0.]

points(S::JacobiSquareSpace,n)=sqrt(fromcanonical(S.domain,gausschebyshev(n,4)[1]))

plan_transform(S::JacobiSquareSpace,n)=gausschebyshev(n,4)

transform(S::JacobiSquareSpace,vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))
function transform(S::JacobiSquareSpace{0,0,0},vals::Vector,xw::(Vector,Vector))
    ## Same as jacobitransform.jl
    x,w=xw
    
    n=length(vals)

    V=jacobip(0:n-1,m+0.5,-0.5,x)'
    nrm=(V.^2)*w
    (V*(w.*vals))./nrm
end

function transform{m}(S::JacobiSquareSpace{m,m,0},vals::Vector,xw::(Vector,Vector))
    ## Same as jacobitransform.jl
    x,w=xw
    
    n=length(vals)

    w2=(1-x).^(m/2)
    mw=w2.*w
    V=jacobip(0:n-int(m/2)-1,m+0.5,-0.5,x)'  
    nrm=(V.^2)*(w2.*mw)    
    (V*(mw.*vals))./nrm
end


#evaluate{m,order}(f::Fun{JacobiSquareSpace{m,order}},x::Number)=x^m*dot(jacobip(0:length(f)-1,m+0.5+order,order-0.5,tocanonical(f,x^2)),f.coefficients)*2^(m/2)



plan_itransform(S::JacobiSquareSpace,n)=points(S,n)
#TODO: general domain
itransform(S::JacobiSquareSpace,cfs::Vector)=itransform(S,cfs,plan_itransform(S,length(cfs)))
itransform{m,a,b}(S::JacobiSquareSpace{m,a,b},cfs::Vector,x)=x.^m.*jacobip(0:length(cfs)-1,a+0.5,b-0.5,tocanonical(S,x.^2))*cfs*2^(m/2)
evaluate{J<:JacobiSquareSpace,T}(f::Fun{J,T},x::Vector)=itransform(f.space,f.coefficients,x)
evaluate{J<:JacobiSquareSpace,T}(f::Fun{J,T},x::Number)=itransform(f.space,f.coefficients,[x])[1]



## Operators



 function Derivative{m,a,b}(S::JacobiSquareSpace{m,a,b})
     # we have D[r^m f(r^2)] = r^{m-1} (m f(r^2) + 2r^2 f'(r^2))
 
     d=domain(S)
     @assert d==Interval(1.,0.)
     
     JS=JacobiSpace(a+.5,b-.5,d)
     D=Derivative(JS)
     M=Multiplication(Fun(identity,d),rangespace(D))
     
     DerivativeWrapper(SpaceOperator((sqrt(2)*m)+(sqrt(2)*2)*M*D,S,JacobiSquareSpace{m-1,a+1,b+1}(d)),1)
end
     

function Derivative(S::JacobiSquareSpace,k::Integer)
     if k==1
         Derivative(S)
     else
         D=Derivative(S)
         DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
     end
end
