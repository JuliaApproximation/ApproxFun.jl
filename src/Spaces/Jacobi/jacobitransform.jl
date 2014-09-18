if isdir(Pkg.dir("FastGaussQuadrature"))
    require("FastGaussQuadrature")
    
    gaussjacobi(n,a,b)=Main.FastGauss.GaussJacobi(n,a,b)
else
    gaussjacobi(n,a,b)=error("Currently require FastGaussQuadrature.jl")    
end

points(S::JacobiSpace,n)=gaussjacobi(n,S.a,S.b)[1]
function transform(S::JacobiSpace,v::Vector)
    n=length(v)
    a=S.a;b=S.b
    x,w=gaussjacobi(n,a,b)
    
    V=jacobip(0:n-1,a,b,x)'
    nrm=(V.^2)*w
    
    V*(w.*v)./nrm
end


evaluate{T}(f::Fun{T,JacobiSpace},x)=dot(jacobip(0:length(f)-1,f.space.a,f.space.b,x),f.coefficients)