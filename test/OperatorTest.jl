using ApproxFun


f=Fun(exp);
d=f.domain;
D=diff(d);
I=eye(d);
Q=integrate(d);

@assert norm((Q+I)*f-(integrate(f)+f)) < eps()
@assert norm((Q)*f-(integrate(f))) < eps()

x=Fun(identity)
X=MultiplicationOperator(x)

@assert norm(ConversionOperator(0:2)\coefficients(x.*f,2)-(x.*f).coefficients) < 100eps()

@assert norm((ConversionOperator(0:2)\(D^2*f))-diff(diff(f))) < 100eps()

@assert norm(X*f-(x.*f)) < 100eps()

@assert norm(ConversionOperator(0:2)*X*f.coefficients-coefficients(x.*f,2)) < 100eps()


## Special functions

x=Fun(identity);
@assert norm(cos(x)-Fun(cos))<10eps()
@assert norm(sin(x)-Fun(sin))<10eps()
@assert norm(exp(x)-Fun(exp))<10eps()
@assert norm(sin(x)./x-Fun(x->sinc(x/π)))<10eps()