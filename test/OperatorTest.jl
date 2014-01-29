using ApproxFun


f=Fun(exp);
d=f.domain;
D=diff(d);
I=eye(d);
Q=integrate(d);

@assert norm((Q+I)*f-(integrate(f)+f)) < eps()
@assert norm((Q)*f-(integrate(f))) < eps()

