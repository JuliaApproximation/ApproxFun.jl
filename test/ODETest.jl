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

u=[Bm;D2-X;Bp]\[airyai(d.a),0.,airyai(d.b)];
@test_approx_eq_eps u(0.) airyai(0.) 10ncoefficients(u)*eps()

u=[D2-X;Bm;Bp]\[0.,airyai(d.a),airyai(d.b)];
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
B=ApproxFun.SpaceOperator(BasisFunctional(floor(w)),Chebyshev(),ApproxFun.ConstantSpace());
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

d=Interval()
D=Derivative(d)
A=D^2-I
κ=nullspace(A)
@test length(κ) == 2

c=[κ(0.);κ'(0.)]\[exp(0.);exp(0.)]
u=(κ*c)[1]
@test_approx_eq u(1.0) e


d=Interval(-50.,5.)
x=Fun(identity,d)
D=Derivative(d)
u=nullspace(D^2-x)
c=[u(d.a); u(d.b)]\[airyai(d.a),airyai(d.b)]
@test norm((u*c)[1]-Fun(airyai,d))<10000eps()


## constant forcing
d = Interval(0.,50.)
D = Derivative(d)
t = Fun(identity,d)

F = D^2 +.5D + I

A= [ 0    ldirichlet(d);
     0    lneumann(d);
     0    rdirichlet(d);
    -1    F; ]

u,x=A\[1.,0.,2.,0.]

@test norm(F*x-u)<1000eps()



## QR tests


S=Chebyshev()
B=dirichlet(S)
D=Derivative(S)

Q,R=qr([B;D^2+I])
u=R\(Q'*[cos(-1.0),cos(1.0)])


@test_approx_eq u(0.) cos(0.0)


S=Chebyshev()
A=[dirichlet(S);Derivative(S)^2 - I]
QR=qrfact(A)
@test_approx_eq (QR\[1.])(0.0) 0.3240271368319427
Q,R=qr(A)
u=(R\(Q'*[1.]))
@test_approx_eq u(0.0)  0.3240271368319427

x=Fun(S)
A=[dirichlet(S);Derivative(S)^2 - exp(im*x)]
QR=qrfact(A)

u=(QR\[1.])
@test_approx_eq u(0.0) (0.3329522068795961 + 0.024616008954634165im)


x=Fun(identity,[-2.,-1.,0.,15.])
sp=space(x)
D=Derivative(sp)
A=[dirichlet(sp);D^2-x]
QR=qrfact(A)
u=QR\[airyai(-2.)]

@test_approx_eq u(0.0) airyai(0.)


## Vector
d=Interval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();

f=Fun(x->[exp(x),cos(x)],d)
A=[B 0;
   Bn 0;
   0 B;
   D^2-I 2.0I;
   0 D+I];


io=ApproxFun.interlace(A)
    co=cache(io)
    n=10

    io=co.op
    ds=domainspace(io)
    rs=rangespace(io)
    di=io.domaininterlacer
    ri=io.rangeinterlacer
    ddims=di.iterator.dimensions
    rdims=ri.iterator.dimensions

    d∞=find(isinf,[ddims...])
    r∞=find(isinf,[rdims...])
    p=length(d∞)


    (l,u)=bandwidths(co.data.bands)
    pad!(co.data,n,n+u)
    co.data
    # r is number of extra rows, ncols is number of extra columns
    r=rank(co.data.fill)
    ncols=mapreduce(d->isfinite(d)?d:0,+,ddims)


    # fill rows
    K=k=1
    while k ≤ r
        if isfinite(dimension(rs[ri[K][1]]))
            co.data.fill.V[co.datasize[2]:end,k] = co.op[K,co.datasize[2]:n+u]
            k += 1
        end
        K += 1
    end

    kr=co.datasize[1]+1:n
    jr=max(1,kr[1]-l):n+u
    io∞=InterlaceOperator(io.ops[r∞,d∞])

    BLAS.axpy!(1.0,view(io∞,kr-r,jr-ncols),view(co.data.bands,kr,jr))
    co.data-io[1:10,1:16]



view(co.op,kr,jr)

bandwidth(view(co.data.bands,kr,jr),2)


kr,jr



kr-r
    jr-ncols

    k=j=1
    njr=(jr[findfirst(ℓ->mod(ℓ-1,p)==j-1,jr-ncols)]:p:jr[end])
    nkr=(kr[findfirst(ℓ->mod(ℓ-1,p)==k-1,kr-r)]:p:kr[end])
    jr1=(njr[1]-ncols+1)÷p
    kr1=(nkr[1]-r+1)÷p
BandedMatrices.banded_axpy!(1.0,view(co.op.ops[r∞[k],d∞[j]],
                    kr1:(kr1+length(nkr)-1),
                    jr1:(jr1+length(njr)-1)),
                    view(co.data.bands,nkr,njr))



nkr
njr


(a,X,S)=(1.0,view(co.op.ops[r∞[k],d∞[j]],
                                        kr1:(kr1+length(nkr)-1),
                                        jr1:(jr1+length(njr)-1)),
                                        view(co.data.bands,nkr,njr))

y
bandwidths(X)
bandwidths(S.parent)
shft=bandshift(S)

bandwidth(view(co.op.ops[r∞[k],d∞[j]],
                    kr1:(kr1+length(nkr)-1),
                    jr1:(jr1+length(njr)-1)),2)


bandwidth(view(co.data.bands,nkr,njr)

view(co.data.bands,nkr,njr)


nkr
njr

kr1
jr1


io.ops[r∞[k],d∞[j]]




nkr


view(co.data.bands,nkr,njr)


for k=1:p,j=1:p
    BLAS.axpy!(1.0,view(co.op.ops[r∞[k],d∞[j]],
                        kr-r,
                        jr-ncols),
    view(co.data.bands,kr,jr))
end

co.datasize=(n,n+u)

QR=qrfact(A)
v=Any[0.,0.,0.,f...]
@test_approx_eq (QR\v)(0.0) [0.0826967758420519,0.5553968826533497]


Q,R=qr(A)
v=Any[0.,0.,0.,f...]
@test_approx_eq (QR\v)(0.0) [0.0826967758420519,0.5553968826533497]
