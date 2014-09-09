using ApproxFun, Base.Test


f=Fun(exp);
d=domain(f);
D=diff(d);
Q=integrate(d);

@test norm((Q+I)*f-(integrate(f)+f)) < eps()
@test norm((Q)*f-(integrate(f))) < eps()

x=Fun(identity)
X=MultiplicationOperator(x)

@test norm(USConversionOperator(0:2,domain(x))\coefficients(x.*f,2)-(x.*f).coefficients) < 100eps()

@test norm((USConversionOperator(0:2,domain(x))\(D^2*f))-diff(diff(f))) < 100eps()

@test norm(X*f-(x.*f)) < 100eps()

@test norm(USConversionOperator(0:2,domain(x))*X*f.coefficients-coefficients(x.*f,2)) < 100eps()


## Special functions

x=Fun(identity);
@test norm(cos(x)-Fun(cos))<10eps()
@test norm(sin(x)-Fun(sin))<10eps()
@test norm(exp(x)-Fun(exp))<10eps()
@test norm(sin(x)./x-Fun(x->sinc(x/π)))<100eps()



## Periodic

a=FFun(t-> cos(t) + sin(2t),3)
d=domain(a)
D=diff(d)
L=D^2+a;
f=FFun(t->cos(cos(t)));
u=L\f;

@test norm(diff(u,2)+a.*u-f) < 10eps()



## Check mixed

d=Interval()
D=diff(d)
x=Fun(identity,d)
A=D*(x*D)
B=D+x*D^2
C=x*D^2+D
@test norm((A-B)[1:10,1:10]|>full)<eps()
@test norm((B-A)[1:10,1:10]|>full)<eps()
@test norm((A-C)[1:10,1:10]|>full)<eps()
@test norm((C-A)[1:10,1:10]|>full)<eps()
@test norm((C-B)[1:10,1:10]|>full)<eps()
@test norm((B-C)[1:10,1:10]|>full)<eps()


