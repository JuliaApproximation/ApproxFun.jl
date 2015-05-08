if isdir(Pkg.dir("FastGaussQuadrature"))
    require("FastGaussQuadrature")
    
    gaussjacobi(n,a,b)=Main.FastGaussQuadrature.gaussjacobi(n,a,b) 
else
    gaussjacobi(n...)=error("Currently require FastGaussQuadrature.jl")    
end



## From FastGaussQuadrature.jl
function gausschebyshev( n::Int64, kind::Int64=1 )
# GAUSS-CHEBYSHEV NODES AND WEIGHTS. 

x = (Array(Float64,n), Array(Float64,n))
# Use known explicit formulas. Complexity O(n).
if kind == 1 
    # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
    x = (cos((2*[n:-1:1]-1)*pi/2n), pi./n*ones(n))
elseif kind == 2 
    # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
    x = (cos([n:-1:1]*pi./(n+1)), pi/(n+1)*sin([n:-1:1]./(n+1)*pi).^2 )
elseif kind == 3 
    # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
    x = (cos(([n:-1:1]-.5)*pi/(n+.5)), 2*pi/(n+.5)*cos(([n:-1:1]-.5)*pi/(2(n+.5))).^2)
elseif kind == 4 
    # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
    x = (cos([n:-1:1]*pi/(n+.5)), 2*pi/(n+.5)*sin([n:-1:1]*pi/(2(n+.5))).^2)
else 
   throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4")) 
end
return x 
end


points(S::Jacobi,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])
function transform(S::Jacobi,v::Vector,xw::@compat(Tuple{Vector,Vector}))
    x,w=xw
    V=jacobip(0:length(v)-1,S.a,S.b,x)'
    nrm=(V.^2)*w
    
    V*(w.*v)./nrm
end

transform(S::Jacobi,v::Vector)=transform(S,v,gaussjacobi(length(v),S.a,S.b))

itransform(S::Jacobi,cfs::Vector,x::Vector)=jacobip(0:length(cfs)-1,S.a,S.b,tocanonical(S,x))*cfs
itransform(S::Jacobi,cfs::Vector)=itransform(S,cfs,points(S,length(cfs)))


evaluate(f::Fun{Jacobi},x::Number)=dot(jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x)),f.coefficients)
evaluate(f::Fun{Jacobi},x::Vector)=jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x))*f.coefficients


## JacobiWeight


function plan_transform(S::JacobiWeight{Jacobi},n::Integer)
    m=S.β
    @assert isapproxinteger(m)
    
    if S.α==S.space.b==0 && S.space.a==2m+1
        gaussjacobi(n,1.,0.)
    elseif S.α==0 && S.space.b ==-0.5 && S.space.a==2m+0.5
        gausschebyshev(n,4) # a=0.5,b==-0.5
    elseif S.α==0 && S.space.b ==-0.5 && S.space.a==2m-0.5
        gausschebyshev(n) # a=-0.5,b==-0.5
    end
end

points(S::JacobiWeight{Jacobi},n)=fromcanonical(S,plan_transform(S,n)[1])


transform(S::JacobiWeight{Jacobi},vals::Vector)=transform(S,vals,plan_transform(S,length(vals)))
function transform(S::JacobiWeight{Jacobi},vals::Vector,xw::@compat(Tuple{Vector,Vector}))
    # Jacobi and JacobiWeight have different a/b orders
    
    if S.α==0 && S.β==0
        return transform(S.space,vals,xw)
    end

    m=S.β
    @assert isapproxinteger(m)
    @assert S.α==0 && (S.space.b==0 && S.space.a==2m+1) || (S.space.b ==-0.5 && S.space.a==2m+0.5)
    
    n=length(vals)
    x,w=xw
    if m==0
        V=jacobip(0:n-1,S.space,x)'
        nrm=(V.^2)*w
        (V*(w.*vals))./nrm
    elseif n>m
        w2=(1-x).^m
        mw=w2.*w
        
        V=jacobip(0:n-int(m)-1,S.space,x)'   
          # only first m coefficients are accurate
          # since Gauss quad is accurate up to polys of degree 2n-1
          # we get one more coefficient because we normalize, so the
          # error for poly of degree 2n is annihilated
        
        
        nrm=(V.^2)*(w2.*mw)
        
        (V*(mw.*vals))./nrm
    else
        [0.]
    end

end



