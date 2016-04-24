using ApproxFun, Base.Test


# This avoids getting killed on Travis.cl
if OS_NAME == :Darwin
    gc_enable(false)
end

## ODEs

d=Interval(-20000.,20000.)
x=Fun(identity,d)
u=[dirichlet(d);Derivative(d)^2+I]\[1.,0.]
u=[dirichlet(d);Derivative(d)^2+I]\[1.,0.]
@time u=[dirichlet(d);Derivative(d)^2+I]\[1.,0.]

println("Cos/Sin: should be ~0.017")


L=Derivative(d)^2-x

@time L[1:12500,1:12500]
@elapsed for k=1:100 L.ops[1][1:1250,1:1250] end
@elapsed for k=1:100 L.ops[2].ops[1][1:1250,1:1250] end
@elapsed for k=1:100 L.ops[2].ops[2][1:1250,1:1250] end


using BandedMatrices
@time (B=bzeros(1250,1250,0,2);)
@time B.data[1,:]=ones(size(B.data,2))
A=L.ops[2].ops[2]
@time for k=1:size(B,1)
    B[k,k] = A[k,k]
end

@time fill!(sub(B.data,1,:),0.5)
B
methods(fill!)

S=sub(A,1:10,1:10)

@which copy(S)

L.ops[2].ops[1].op.ops[2]

@elapsed for k=1:100 L.ops[2].ops[1].op.ops[1][1:1250,1:1250] end
@elapsed for k=1:100 L.ops[2].ops[1].op.ops[2][1:1250,1:1250] end

@which L.ops[2].ops[1].op.ops[2][1,3]


?Diagonal


d=Interval(-1000.,5.)
x=Fun(identity,d)
u=[dirichlet(d);Derivative(d)^2-x]\[1.,0.]
u=[dirichlet(d);Derivative(d)^2-x]\[1.,0.]
@time u=[dirichlet(d);Derivative(d)^2-x]\[1.,0.]
println("Airy: should be ~0.013")
u(0.)
@profile u=[dirichlet(d);Derivative(d)^2-x]\[1.,0.]
Profile.print()
S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+(7+2x+6x^2)
B=dirichlet(S)
n=20000
rhs=ones(n+2)
u=[B;L]\rhs
u=[B;L]\rhs
@time u=[B;L]\rhs
println("Poly: should be ~0.025")


S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+cos(x)
B=dirichlet(S)
n=2000
rhs=ones(n+2)
u=linsolve([B;L],rhs;maxlength=Inf)
u=linsolve([B;L],rhs;maxlength=Inf)
@time u=linsolve([B;L],rhs;maxlength=Inf)
println("Cos: should be ~0.0075")

S=Chebyshev()
x=Fun(identity,S)
D=Derivative(S)
L=D^2+sin(x)
B=dirichlet(S)
n=2000
rhs=ones(n+2)
u=linsolve([B;L],rhs;maxlength=Inf)
u=linsolve([B;L],rhs;maxlength=Inf)
@time u=linsolve([B;L],rhs;maxlength=Inf)
println("Sin: should be ~0.011")


## Piecewise

x=Fun(identity,[-20.,-10.,-5.,0.,1.,15.])
sp=space(x)
D=Derivative(sp)

u=[dirichlet(sp);
    D^2-x]\[airyai(-10.)];
@time u=[dirichlet(sp);
    D^2-x]\[airyai(-10.)];

println("Piecewise Airy: should be ~0.008")


## Vector
d=Interval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();

f=Fun(x->[exp(x),cos(x)],d)

A=[B 0;
   Bn 0;
   0 B;
   D^2-I 2.I;
   0 D+I];

u=A\Any[0.,0.,0.,f]
@time u=A\Any[0.,0.,0.,f]
println("Systems: should be ~0.0008")


d=Interval(-300.,5.)
x=Fun(identity,d)
A=Derivative(d)^2-x
u=null(A)
u=null(A)
@time u=null(A)
println("Null Airy: should be ~0.061")


if OS_NAME == :Darwin
    gc_enable(true)
end
