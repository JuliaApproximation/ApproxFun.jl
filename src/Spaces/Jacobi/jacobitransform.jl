if isdir(Pkg.dir("FastGaussQuadrature"))
    require("FastGaussQuadrature")
    
    gaussjacobi(n,a,b)=Main.FastGaussQuadrature.GaussJacobi(n,a,b)
else
    gaussjacobi(n,a,b)=error("Currently require FastGaussQuadrature.jl")    
end

points(S::JacobiSpace,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])
function transform(S::JacobiSpace,v::Vector)
    n=length(v)
    a=S.a;b=S.b
    x,w=gaussjacobi(n,a,b)
    
    V=jacobip(0:n-1,a,b,x)'
    nrm=(V.^2)*w
    
    V*(w.*v)./nrm
end

function itransform(S::JacobiSpace,cfs::Vector)
    n=length(cfs)
    x=points(S,n)
    jacobip(0:n-1,S.a,S.b,x)*cfs
end


evaluate(f::Fun{JacobiSpace},x)=dot(jacobip(0:length(f)-1,f.space.a,f.space.b,x),f.coefficients)


## JacobiWeightSpace

function points(S::JacobiWeightSpace{JacobiSpace},n)
    m=S.β
    if S.α==S.space.b==0 && S.space.a==2m+1
        gaussjacobi(n,1.,0.)[1]
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

