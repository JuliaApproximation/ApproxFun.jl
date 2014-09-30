if isdir(Pkg.dir("FastGaussQuadrature"))
    require("FastGaussQuadrature")
    
    gaussjacobi(n,a,b)=Main.FastGaussQuadrature.gaussjacobi(n,a,b)
else
    gaussjacobi(n,a,b)=error("Currently require FastGaussQuadrature.jl")    
end

points(S::JacobiSpace,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])
function transform(S::JacobiSpace,v::Vector,x::Vector,w::Vector)
    V=jacobip(0:length(v)-1,S.a,S.b,x)'
    nrm=(V.^2)*w
    
    V*(w.*v)./nrm
end

transform(S::JacobiSpace,v::Vector)=transform(S,v,gaussjacobi(length(v),S.a,S.b)...)

itransform(S::JacobiSpace,cfs::Vector,x::Vector)=jacobip(0:length(cfs)-1,S.a,S.b,tocanonical(S,x))*cfs
itransform(S::JacobiSpace,cfs::Vector)=itransform(S,cfs,points(JacobiSpace(S.a,S.b),length(cfs)))


evaluate(f::Fun{JacobiSpace},x)=dot(jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x)),f.coefficients)


## JacobiWeightSpace

function points(S::JacobiWeightSpace{JacobiSpace},n)
    m=S.β
    if S.α==S.space.b==0 && S.space.a==2m+1
        fromcanonical(S,gaussjacobi(n,1.,0.)[1])
    else
        error("JacobiWeightSpace{JacobiSpace} only implemented for special case a=2m+1,b=0 currently")
    end
end
function transform(S::JacobiWeightSpace{JacobiSpace},vals::Vector)
    # JacobiSpace and JacobiWeightSpace have different a/b orders
    m=S.β
    if S.α==S.space.b==0 && S.space.a==2m+1
        n=length(vals)
        x,w=gaussjacobi(n,1.,0.)
        w2=(1-x).^m
        mw=w2.*w
        
        V=jacobip(0:n-1,S.space,x)'
        
        
        nrm=(V.^2)*(w2.*mw)
        
        (V*(mw.*vals))./nrm
    else
        error("JacobiWeightSpace{JacobiSpace} only implemented for special case a=2m+1,b=0 currently")
    end    
end

