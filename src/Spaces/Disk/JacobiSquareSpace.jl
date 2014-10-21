
# represents function of the form r^m g(r^2)
# as r^m P^{m+.5+order,m-.5+order}(r^2)

immutable JacobiSquareSpace{m,order} <: IntervalDomainSpace
    domain::Union(Interval,AnyDomain)
end
JacobiSquareSpace(m,order,d)=JacobiSquareSpace{m,order}(d)
JacobiSquareSpace(m,d::Domain)=JacobiSquareSpace(m,0,d)
JacobiSquareSpace(m::Integer,o::Integer)=JacobiSquareSpace(m,o,Interval(1.,0.))
JacobiSquareSpace(m::Integer)=JacobiSquareSpace(m,0)
spacescompatible{m,order}(a::JacobiSquareSpace{m,order},b::JacobiSquareSpace{m,order})=true

# We assume domains is [1.,0.]

points(S::JacobiSquareSpace,n)=sqrt(fromcanonical(S.domain,gausschebyshev(n,4)[1]))

plan_transform(S::JacobiSquareSpace,n)=gausschebyshev(n,4)

transform{m}(S::JacobiSquareSpace{m,0},vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))
function transform(S::JacobiSquareSpace{0,0},vals::Vector,xw::(Vector,Vector))
    ## Same as jacobitransform.jl
    x,w=xw
    
    n=length(vals)

    V=jacobip(0:n-1,m+0.5,-0.5,x)'
    nrm=(V.^2)*w
    (V*(w.*vals))./nrm
end

function transform{m}(S::JacobiSquareSpace{m,0},vals::Vector,xw::(Vector,Vector))
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
itransform{m,order}(S::JacobiSquareSpace{m,order},cfs::Vector,x)=x.^m.*jacobip(0:length(cfs)-1,m+2order+0.5,order-0.5,tocanonical(S,x.^2))*cfs*2^(m/2)
evaluate{J<:JacobiSquareSpace,T}(f::Fun{J,T},x::Vector)=itransform(f.space,f.coefficients,x)
evaluate{J<:JacobiSquareSpace,T}(f::Fun{J,T},x::Number)=itransform(f.space,f.coefficients,[x])[1]



## Operators



 function Derivative{m,order}(S::JacobiSquareSpace{m,order})
     # we have D[r^m f(r^2)] = r^{m-1} (m f(r^2) + 2r^2 f'(r^2))
 
     d=domain(S)
     @assert d==Interval(1.,0.)
     
     JS=JacobiSpace(m+2order+.5,order-.5,d)
     D=Derivative(JS)
     M=Multiplication(Fun(identity,d),rangespace(D))
     
     DerivativeWrapper(SpaceOperator((sqrt(2)*m)+(sqrt(2)*2)*M*D,S,JacobiSquareSpace{m-1,order+1}(d)),1)
end
     
     
#     
#     r=Fun(identity,d)
#     
#     #multiplying by r^2 is same as multiplying Jacobi series by r
#     
#     
#     
#     
# 
#     if S.α==S.β==0
#         DerivativeWrapper(SpaceOperator(Derivative(S.space),S,JacobiWeightSpace(0.,0.,rangespace(Derivative(S.space)))),1)
#     elseif S.α==0
#         x=Fun(identity,d)
#         M=tocanonical(d,x)
#         Mp=tocanonicalD(d,d.a)            
#         DD=(-Mp*S.β)*I +(1-M)*Derivative(S.space)
#         DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(0.,S.β-1,rangespace(DD))),1)
#     elseif S.β==0
#         x=Fun(identity,d)
#         M=tocanonical(d,x)
#         Mp=tocanonicalD(d,d.a)        
#         DD=(Mp*S.α)*I +(1+M)*Derivative(S.space)
#         DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(S.α-1,0.,rangespace(DD))),1)
#     else 
#         x=Fun(identity,d)
#         M=tocanonical(d,x)
#         Mp=tocanonicalD(d,d.a)
#         DD=(Mp*S.α)*(1-M) - (Mp*S.β)*(1+M) +(1-M.^2)*Derivative(S.space)
#         DerivativeWrapper(SpaceOperator(DD,S,JacobiWeightSpace(S.α-1,S.β-1,rangespace(DD))),1)
#     end
# 
# end
# 
# function Derivative(S::JacobiWeightSpace,k::Integer)
#     if k==1
#         Derivative(S)
#     else
#         D=Derivative(S)
#         DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
#     end
# end
