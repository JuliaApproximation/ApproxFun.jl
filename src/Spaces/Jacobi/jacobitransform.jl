if isdir(Pkg.dir("FastGaussQuadrature"))
    import FastGaussQuadrature

    gaussjacobi(n,a,b)=Main.FastGaussQuadrature.gaussjacobi(n,a,b)
else
    ## From FastGaussQuarture.jl
    # Only the ±0.5,±0.5 are used
    function gaussjacobi(n::Int64, a, b)
    #GAUSS-JACOBI QUADRATURE NODES AND WEIGHTS

        if ( a == 0 && b == 0 )
            x = gausslegendre(n)
        elseif ( a == -0.5 && b == -0.5 )
            x = gausschebyshev(n,1)
        elseif ( a == 0.5 && b == 0.5)
            x = gausschebyshev(n,2)
        elseif ( a == -0.5 && b == 0.5)
            x = gausschebyshev(n,3)
        elseif ( a == 0.5 && b == -0.5 )
            x = gausschebyshev(n,4)
        elseif ( n == 0 )
            x = (Float64[],Float64[])
        elseif ( n == 1 )
            x = ([(b-a)/(a+b+2)], [2^(a+b+1)*beta(a+1, b+1)])
        elseif ( n <= 100 )
            x = JacobiRec(n, a, b)
        elseif ( n > 100 )
            x = JacobiAsy(n, a, b)
        else
            error("1st argument must be a positive integer.")
        end
        return x
    end

end



## From FastGaussQuadrature.jl
function gausschebyshev( n::Int64, kind::Int64=1 )
    # GAUSS-CHEBYSHEV NODES AND WEIGHTS.

    # Use known explicit formulas. Complexity O(n).
    if kind == 1
        # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        ([cos((2*k-1)*pi/2n) for k=n:-1:1], fill(pi./n,n))
    elseif kind == 2
        # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
        ([cos(k*pi./(n+1)) for k=n:-1:1], [pi/(n+1)*sin(k./(n+1)*pi).^2 for k=n:-1:1])
    elseif kind == 3
        # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        ([cos((k-.5)*pi/(n+.5)) for k=n:-1:1], [2*pi/(n+.5)*cos((k-.5)*pi/(2(n+.5))).^2 for k=n:-1:1])
    elseif kind == 4
        # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        ([cos(k*pi/(n+.5)) for k=n:-1:1], [2*pi/(n+.5)*sin(k*pi/(2(n+.5))).^2 for k=n:-1:1])
    else
       throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4"))
    end
end


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

evaluate(f::Fun{Jacobi},x::Number)=length(f)==0?zero(x):dot(jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x)),f.coefficients)
evaluate(f::Fun{Jacobi},x::Vector)=jacobip(0:length(f)-1,f.space.a,f.space.b,tocanonical(f,x))*f.coefficients


## JacobiWeight


plan_transform(S::JacobiWeight{Jacobi},v::Vector)=gaussjacobi(length(v),S.β,S.α)

points(S::JacobiWeight{Jacobi},n)=fromcanonical(S,plan_transform(S,ones(n))[1])

function transform(S::JacobiWeight{Jacobi},vals::Vector,plan::@compat(Tuple{Vector,Vector}))
    # Jacobi and JacobiWeight have different a/b orders

    if S.α==0 && S.β==0
        return transform(S.space,vals,plan)
    end

    m=S.β
    @assert isapproxinteger(m)
    @assert S.α==0 && (S.space.b==0 && S.space.a==2m+1) || (S.space.b ==-0.5 && S.space.a==2m+0.5)

    n=length(vals)
    x,w=plan
    if m==0
        V=jacobip(0:n-1,S.space,x)'
        nrm=(V.^2)*w
        (V*(w.*vals))./nrm
    elseif n>m
        w2=(1-x).^m
        mw=w2.*w

        V=jacobip(0:n-round(Int,m)-1,S.space,x)'
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
