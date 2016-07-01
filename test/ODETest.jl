using ApproxFun,Base.Test
import ApproxFun.Multiplication


##Airy equation


d=Interval(-10.,5.);
S=Chebyshev(d)


Bm=Evaluation(d,d.a);
Bp=Evaluation(d,d.b);
B=[Bm;Bp];
D2=Derivative(d,2);
X=Multiplication(Fun(x->x,d));

u=[B;D2-X]\[airyai(d.a),airyai(d.b),0.];

@test_approx_eq_eps u(0.) airyai(0.) 10ncoefficients(u)*eps()



d=Interval(-1000.,5.);
Bm=Evaluation(d,d.a);
Bp=Evaluation(d,d.b);
B=[Bm;Bp];
D2=Derivative(d,2);
X=Multiplication(Fun(x->x,d));

u=[B;D2-X]\[airyai(d.a),airyai(d.b),0.];
@test_approx_eq_eps u(0.) airyai(0.) 10ncoefficients(u)*eps()



B=neumann(d);
A=[B;D2-X];
b=[airyaiprime(d.a),airyaiprime(d.b),0.];

u=A\b;

@test_approx_eq_eps u(0.) airyai(0.) 10ncoefficients(u)*eps()

## Neumann condition



f=Fun(x->x.^2)
D=Derivative(domain(f))
@test norm(D*f-f')<100eps()


##Test versus exp

f=Fun(x->-x.^2)
g=Fun(t->exp(-t.^2))

@test norm(Fun(t->exp(f(t)))-g)<= 100eps()

fp=f';
Bm=Evaluation(domain(f),domain(f).a);
u=[Bm,Derivative(domain(f)) - fp]\[exp(f(domain(f).a)),0.];
@test norm(u-g)<100eps()



## Oscillatory integral

f=Fun(exp);
D=Derivative(domain(f));
w=10.;
B=BasisFunctional(floor(w));
A=[B;D+1im*w*I];
u = A\[0.,f];
@test_approx_eq (u(1.)exp(1im*w)-u(-1.)exp(-1im*w)) (-0.18575766879136255 + 0.17863980562549928im )



## Bessel

d=Interval()
D=Derivative(d)
x=Fun(identity,d)
A=x^2*D^2+x*D+x^2
u=[dirichlet(d)[1];A]\[besselj(0,d.a),0.];

@test_approx_eq u(0.1) besselj(0.,0.1)
@test norm(A*u)<10eps()
@test norm(Fun(A.ops[1]*u,d)-x.^2.*differentiate(u,2))<eps()
@test norm(Fun(A.ops[2]*u,d)-x.*u') < eps()
@test norm(Fun(A.ops[end]*u,d)-x.^2.*u) < eps()
@test norm(x.^2.*u'' + x.*u' + x.^2.*u)<10eps()





## Null space

d=Interval(-50.,5.)
D=Derivative(d)
x=Fun(identity,d)
u=nullspace(D^2-x)
c=[evaluate(u,d.a)'; evaluate(u,d.b)']\[airyai(d.a),airyai(d.b)]
@test norm(dot(c,u)-Fun(airyai,d))<eps(1000.)






## constant forcing


d = Interval(0.,50.)
D = Derivative(d)
t = Fun(identity,d)

F = D^2 +.5D + 1

A= [0  ldirichlet(d);
    0    lneumann(d);
    0    rdirichlet(d);
    -1    F; ]

u,x=A\[1.,0.,2.,0.]

@test norm(F*x-u)<1000eps()
