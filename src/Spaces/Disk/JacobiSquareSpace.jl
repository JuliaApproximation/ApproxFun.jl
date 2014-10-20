
# represents function of the form r^m g(r^2)


immutable JacobiSquareSpace <: IntervalDomainSpace
    m::Int
    domain::Union(Interval,AnyDomain)
end

JacobiSquareSpace(m)=JacobiSquareSpace(m,Interval(1.,0.))
spacescompatible(a::JacobiSquareSpace,b::JacobiSquareSpace)=a.m==b.m

# We assume domains is [1.,0.]

points(S::JacobiSquareSpace,n)=sqrt(fromcanonical(S.domain,gausschebyshev(n,4)[1]))

plan_transform(S::JacobiSquareSpace,n)=gausschebyshev(n,4)

transform(S::JacobiSquareSpace,vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))
function transform(S::JacobiSquareSpace,vals::Vector,xw::(Vector,Vector))
    ## Same as jacobitransform.jl
    x,w=xw
    m=S.m
    n=length(vals)
    if m==0
        V=jacobip(0:n-1,m+0.5,-0.5,x)'
        nrm=(V.^2)*w
        (V*(w.*vals))./nrm
    else    
        w2=(1-x).^(m/2)
        mw=w2.*w
        V=jacobip(0:n-int(m/2)-1,m+0.5,-0.5,x)'  
        nrm=(V.^2)*(w2.*mw)    
        (V*(mw.*vals))./nrm
    end
end

#TODO: general domain
evaluate(f::Fun{JacobiSquareSpace},x::Number)=x^f.space.m*dot(jacobip(0:length(f)-1,f.space.m+0.5,-0.5,tocanonical(f,x^2)),f.coefficients)*2^(f.space.m/2)



plan_itransform(S::JacobiSquareSpace,n)=points(S,n)
itransform(S::JacobiSquareSpace,cfs::Vector)=itransform(S,cfs,plan_itransform(S,length(cfs)))
itransform(S::JacobiSquareSpace,cfs::Vector,x)=x.^S.m.*jacobip(0:length(cfs)-1,S.m+0.5,-0.5,tocanonical(S,x.^2))*cfs*2^(S.m/2)
evaluate(f::Fun{JacobiSquareSpace},x::Vector)=itransform(f.space,f.coefficients,x)