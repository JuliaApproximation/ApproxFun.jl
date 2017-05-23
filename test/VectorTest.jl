using ApproxFun,Base.Test
    import ApproxFun:interlace,Multiplication,ConstantSpace,PointSpace,ArraySpace,testblockbandedoperator


## Vector*Vector{Fun}

f = Fun(x->[exp(x),cos(x)])
@test ([1 2]*f)(0.1) ≈ [1 2]*f(0.1)
@test_broken f.'*[1,2] ≈ f(0.1).'*[1,2]

## Matrix*Vector{Fun}

a = [1 2; 3 4]
# Chebyshev Vector
@test (a*f)(0.1) ≈ [exp(0.1)+2cos(0.1); 3exp(0.1)+4cos(0.1)]
@test (a*f)(0.1) ≈ a*f(0.1)
@test Fun(a)*f ≈ a*f
@test Fun(a*Array(f)) ≈ a*f

@test (f+1)(0.1) ≈ f(0.1)+1


# Chebyshev Matrix
m = Fun(x->[exp(x) cos(x); sin(x) airyai(x)])
@test (a*m)(0.1) ≈ [exp(0.1)+2sin(0.1) cos(0.1)+2airyai(0.1);
                    3exp(0.1)+4sin(0.1) 3cos(0.1)+4airyai(0.1)]
@test (a*m)(0.1) ≈ a*m(0.1)
@test (m*a)(0.1) ≈ m(0.1)*a
@test Fun(a)*m   ≈ a*m
@test m*Fun(a)   ≈ m*a
@test Fun(a*Array(m))   ≈ a*m
@test Fun(Array(m)*a)   ≈ m*a

@test (a+m)(0.1) ≈ a+m(0.1)
@test (m+a)(0.1) ≈ m(0.1)+a

@test (m+1)(0.1) ≈ m(0.1)+1
@test (m+I)(0.1) ≈ m(0.1)+I


# CosSpace Vector
f = Fun(θ->[1,cos(θ)],CosSpace())
@test (a*f)(0.1) ≈ [1+2cos(0.1); 3+4cos(0.1)]
@test (a*f)(0.1) ≈ a*f(0.1)
@test Fun(a)*f ≈ a*f
@test Fun(a*Array(f)) ≈ a*f

@test (f+1)(0.1) ≈ f(0.1)+1

# CosSpace Matrix
m = Fun(θ->[1 cos(θ); cos(2θ) cos(cos(θ))],CosSpace())
@test (a*m)(0.1) ≈ a*m(0.1)
@test (m*a)(0.1) ≈ m(0.1)*a
@test Fun(a)*m   ≈ a*m
@test Fun(a*Array(m))   ≈ a*m

@test (a+m)(0.1) ≈ a+m(0.1)
@test (m+a)(0.1) ≈ m(0.1)+a

@test (m+1)(0.1) ≈ m(0.1)+1
@test (m+I)(0.1) ≈ m(0.1)+I


# SinSpace Vector
f = Fun(θ->[sin(θ),sin(2θ)],SinSpace())
@test (a*f)(0.1) ≈ a*f(0.1)
@test Fun(a)*f ≈ a*f
@test Fun(a*Array(f)) ≈ a*f

@test all(sp -> sp isa SinSpace, space(a*f).spaces)

@test (f+1)(0.1) ≈ f(0.1)+1

# SinSpace Matrix
# CosSpace Matrix
m = Fun(θ->[sin(3θ) sin(θ); sin(2θ) sin(sin(θ))],SinSpace())
@test (a*m)(0.1) ≈ a*m(0.1)
@test (m*a)(0.1) ≈ m(0.1)*a
@test Fun(a)*m   ≈ a*m
@test Fun(a*Array(m))   ≈ a*m

@test all(sp -> sp isa SinSpace, space(a*m).spaces)

@test (a+m)(0.1) ≈ a+m(0.1)
@test (m+a)(0.1) ≈ m(0.1)+a

@test (m+1)(0.1) ≈ m(0.1)+1
@test (m+I)(0.1) ≈ m(0.1)+I

# Fourier Vector
f = Fun(θ->[sin(θ),sin(2θ)],Fourier())
@test (a*f)(0.1) ≈ a*f(0.1)
@test Fun(a)*f ≈ a*f
@test Fun(a*Array(f)) ≈ a*f




## Matrix{Fun}*Matrix{Fun}

x=Fun()

A = [x x; x x]

@test norm(map(norm,A*A-[2x^2 2x^2; 2x^2 2x^2])) <eps()

@test norm((A*Fun(A)-Fun([2x^2 2x^2; 2x^2 2x^2])).coefficients) < eps()
@test norm((Fun(A)*Fun(A)-Fun([2x^2 2x^2; 2x^2 2x^2])).coefficients) < eps()





d=Interval()
D=Derivative(d);
B=ldirichlet();
Bn=lneumann();
A=[B 0;
   0 B;
   D-I 2I;
   0I D+I];

f=Fun(x->[exp(x),cos(x)],d)


b=Any[0.,0.,f...]
f1,f2=vec(f)
@time u=A\b
u1=vec(u)[1];u2=vec(u)[2];

@test norm(u1'-u1+2u2-f1)<10eps()
@test norm(u2'+u2-f2)<10eps()

Ai=Operator(A)
@time u=Ai\b
u1=vec(u)[1];u2=vec(u)[2];

@test norm(u1'-u1+2u2-f1)<10eps()
@test norm(u2'+u2-f2)<10eps()





A=[B 0;
   Bn 0;
   0 B;
   D^2-I 2I;
   0 D+I];


b=Any[0.,0.,0.,f...]


@time u=A\b
u1=vec(u)[1];u2=vec(u)[2];


@test norm(differentiate(u1,2)-u1+2u2-f1)<2eps()
@test norm(u2'+u2-f2)<2eps()

Ai=Operator(A)
@time u=Ai\b
u1=vec(u)[1];u2=vec(u)[2];


@test norm(u1''-u1+2u2-f1)<2eps()
@test norm(u2'+u2-f2)<2eps()






## Matrix exponential

n=4
d=fill(Interval(0.,1.),n)
B=Evaluation(d,0.)
D=Derivative(d)
A=rand(n,n)
L=[B;D-A]
@time u=L\eye(2n,n)
@test norm(u(1.)-expm(A))<eps(1000.)

n=4
d=fill(Interval(0.,1.),n)
B=Evaluation(d,0.)
D=Derivative(d)
A=rand(n,n)
L=[B;D-A]
@time u=L\eye(2n,2)
@test norm(u(1.)-expm(A)[:,1:2])<eps(1000.)




## Multiplication

d = Interval()
t=Fun(identity,d)
f = Fun([t^2, sin(t)])
@test norm(((Derivative(space(f))*f)-Fun(t->[2t,cos(t)])).coefficients)<100eps()
@test norm((([1 2;3 4]*f)-Fun(t->[t^2+2sin(t),3t^2+4sin(t)])).coefficients)<100eps()



## Multiplication operator

Γ=Circle() ∪ Circle(0.5)


f=Fun(z->in(z,component(Γ,2)) ? 1:z,Γ)
@test f(exp(0.1im)) ≈ exp(0.1im)
@test f(0.5exp(0.1im)) ≈ 1


G=Fun(z->in(z,component(Γ,2))?[1 -z^(-1); 0 1]:
                   [z 0; 0 z^(-1)],Γ);


@test G(exp(0.1im)) ≈ [exp(0.1im) 0 ; 0 exp(-0.1im)]
@test G(0.5exp(0.1im)) ≈ [1 -2exp(-0.1im) ; 0 1]



G1=Fun(Array(G)[:,1])

@test G1(exp(0.1im)) ≈ [exp(0.1im),0.]
@test G1(0.5exp(0.1im)) ≈ [1,0.]

M=Multiplication(G,space(G1))
testblockbandedoperator(M)


for z in (0.5exp(0.1im),exp(0.2im))
    @test G[1,1](z) ≈ G[1](z)
    @test (M.op.ops[1,1]*G1[1])(z) ≈ M.f[1,1](z)*G1[1](z)
    @test (M.op.ops[2,1]*G1[1])(z) ≈ M.f[2,1](z)*G1[1](z)
    @test (M.op.ops[1,2]*G1[2])(z) ≈ M.f[1,2](z)*G1[2](z)
    @test (M.op.ops[2,2]*G1[2])(z) ≈ M.f[2,2](z)*G1[2](z)
end


u=M*G1
@test norm(u(exp(.1im))-[exp(.2im),0])<100eps()
@test norm(u(.5exp(.1im))-[1,0])<100eps()


# Vector operations
@test (Fun(x->[1., 2.]) + [2, 2])(0.) ≈ [3., 4.]



G=Fun(z->[-1 -3; -3 -1]/z +
         [ 2  2;  1 -3] +
         [ 2 -1;  1  2]*z,Circle())


@test G[1,1](exp(0.1im)) == G(exp(0.1im))[1,1]


F̃ = (G-I)[:,1]
F=Fun((G-I)[:,1])

@test F(exp(0.1im)) ≈ [-exp(-0.1im)+1+2exp(0.1im);-3exp(-0.1im)+1+1exp(0.1im)]
@test Fun(F̃,space(F))(exp(0.1im)) ≈ [-exp(-0.1im)+1+2exp(0.1im);-3exp(-0.1im)+1+1exp(0.1im)]

@test coefficients(F̃,space(F)) == F.coefficients
@test Fun(F̃,space(F)) == F

@test F==Fun(vec(F),space(F))

@test inv(G(exp(0.1im))) ≈ inv(G)(exp(0.1im))
@test Fun(eye(2),space(G))(exp(0.1im)) ≈ eye(2)
@test Fun(I,space(G))(exp(0.1im)) ≈ eye(2)


## Check conversion

f=Fun(t->[cos(t) 0;sin(t) 1],-π..π)
g=Fun(f,Space(PeriodicInterval(-π,π)))
@test g(.1) ≈ f(.1)



## Interlace test
S1=Chebyshev()^2
S2=Chebyshev()
TS=ArraySpace([ConstantSpace(),S1,ConstantSpace(),S2,PointSpace([1.,2.])])
f=Fun(TS,collect(1:10))
@test f[1] == Fun(TS[1],[1.])
@test f[2] == Fun(TS[2],[2.,6.,7.,10.])
@test f[3] == Fun(TS[3],[3.])
@test f[4] == Fun(TS[4],[4.,8.])
@test f[5] == Fun(TS[5],[5.,9.])

## Operator * Matrix

D=Derivative()

u=D*[Fun(exp) Fun(cos)]
@test u(0.1) ≈ [exp(0.1) -sin(0.1)]
