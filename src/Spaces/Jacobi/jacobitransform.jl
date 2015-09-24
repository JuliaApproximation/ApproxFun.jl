
points(S::Jacobi,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])

plan_transform(S::Jacobi,v::Vector) = gaussjacobi(length(v),S.a,S.b)
plan_itransform(S::Jacobi,cfs::Vector) = points(S,length(cfs))
function transform(S::Jacobi,vals,plan::@compat(Tuple{Vector,Vector}))
    x,w = plan
    V=jacobip(0:length(vals)-1,S.a,S.b,x)'
    nrm=(V.^2)*w

    V*(w.*vals)./nrm
end
itransform(S::Jacobi,cfs,plan::Vector) = jacobip(0:length(cfs)-1,S.a,S.b,tocanonical(S,plan))*cfs

evaluate{J<:Jacobi}(f::Fun{J},x::Number)=length(f)==0?zero(x):dot(jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x)),f.coefficients)
evaluate{J<:Jacobi}(f::Fun{J},x::Vector)=jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x))*f.coefficients
