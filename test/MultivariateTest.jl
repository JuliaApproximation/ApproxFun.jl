using ApproxFun, Base.Test
    import ApproxFun: testbandedblockbandedoperator


@time for k=0:5,j=0:5
    ff=(x,y)->cos(k*acos(x))*cos(j*acos(y))
    f=Fun(ff,Interval()^2)
    @test f(0.1,0.2) ≈ ff(0.1,0.2)
end

@time for k=0:5,j=0:5
    ff=(x,y)->cos(k*acos(x/2))*cos(j*acos(y/2))
    f=Fun(ff,Interval(-2,2)^2)
    @test f(0.1,0.2) ≈ ff(0.1,0.2)
end


@time for k=0:5,j=0:5
    ff=(x,y)->cos(k*acos(x-1))*cos(j*acos(y-1))
    f=Fun(ff,Interval(0,2)^2)
    @test f(0.1,0.2) ≈ ff(0.1,0.2)
end


## Try constructor variants

ff=(x,y)->exp(-10(x+.2)^2-20(y-.1)^2)*cos(x*y)
gg=x->exp(-10(x[1]+.2)^2-20(x[1]-.1)^2)*cos(x[1]*x[2])
f=Fun(ff,Interval()^2,10000)
@test f(0.,0.) ≈ ff(0.,0.)

f=Fun(gg,Interval()^2,10000)
@test f(0.,0.) ≈ ff(0.,0.)

f=Fun(ff,Interval()^2)
@test f(0.,0.) ≈ ff(0.,0.)
f=Fun(gg,Interval()^2)
@test f(0.,0.) ≈ ff(0.,0.)


f=Fun(ff)
@test f(0.,0.) ≈ ff(0.,0.)
f=Fun(gg)
@test f(0.,0.) ≈ ff(0.,0.)


# Fun +-* constant
f=Fun((x,y)->exp(x)*cos(y))

@test f(0.1,0.2)+2 ≈ (f+2)(0.1,0.2)
@test f(0.1,0.2)-2 ≈ (f-2)(0.1,0.2)
@test f(0.1,0.2)*2 ≈ (f*2)(0.1,0.2)


## LowRankFun

@time F = LowRankFun((x,y)->besselj0(10(y-x)),Chebyshev(),Chebyshev())

@test F(.123,.456) ≈ besselj0(10(.456-.123))

@time G = LowRankFun((x,y)->besselj0(10(y-x));method=:Cholesky)

@test G(.357,.246) ≈ besselj0(10(.246-.357))





## 1D in 2D

d=Segment((0.,0.),(1.,1.))
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d)
@test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))
@test f(ApproxFun.Vec(0.5,0.5)) ≈ exp(-0.5-2cos(0.5))

f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d,20)
@test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d)
@test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d,20)
@test f(0.5,0.5) ≈ exp(-0.5-2cos(0.5))


d=Circle((0.,0.),1.)
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),Fourier(d),40)
@test f(cos(0.1),sin(0.1)) ≈ exp(-cos(0.1)-2cos(sin(0.1)))
@test f(ApproxFun.Vec(cos(0.1),sin(0.1))) ≈ exp(-cos(0.1)-2cos(sin(0.1)))

f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d),40)
@test f(cos(0.1),sin(0.1)) ≈ exp(-cos(0.1)-2cos(sin(0.1)))


f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d))
@test f(cos(0.1),sin(0.1)) ≈ exp(-cos(0.1)-2cos(sin(0.1)))

println("    Calculus tests")

## Sum

ff=(x,y)->(x-y)^2*exp(-x^2/2.-y^2/2)
f=Fun(ff,Domain(-4..4)^2)
@test f(0.1,0.2) ≈ ff(0.1,0.2)

@test sum(f,1)(0.1) ≈ 2.5162377980828357
f=LowRankFun(f)
@test evaluate(f.A,0.1) ≈ map(f->f(0.1),f.A)


## Kron operator
Mx=Multiplication(Fun(cos),Chebyshev())
My=Multiplication(Fun(sin),Chebyshev())
K=Mx⊗My

@test ApproxFun.BandedBlockBandedMatrix(view(K,1:10,1:10)) ≈ [K[k,j] for k=1:10,j=1:10]
C=Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical(1)⊗Ultraspherical(1))
@test C[1:100,1:100] ≈ Float64[C[k,j] for k=1:100,j=1:100]


@time let d = Space(0..1) * Space(0..2)
    Dx = Derivative(d, [1,0])
    f = Fun((x,y) -> sin(x) * cos(y), d)
    fx = Fun((x,y) -> cos(x) * cos(y), d)
    @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
    Dy = Derivative(d, [0,1])
    fy = Fun((x,y) -> -sin(x) * sin(y), d)
    @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
    L=Dx+Dy
    @test (L*f)(0.2,0.3) ≈ (fx(0.2,0.3)+fy(0.2,0.3))

    B=ldirichlet(d[1])⊗ldirichlet(d[2])
    @test abs(Number(B*f)-f(0.,0.)) ≤ 10eps()

    B=Evaluation(d[1],0.1)⊗ldirichlet(d[2])
    @test Number(B*f) ≈ f(0.1,0.)

    B=Evaluation(d[1],0.1)⊗Evaluation(d[2],0.3)
    @test Number(B*f) ≈ f(0.1,0.3)

    B=Evaluation(d,(0.1,0.3))
    @test Number(B*f) ≈ f(0.1,0.3)
end



## x,y constructor

@time let d=Interval()^2
    x,y=Fun(d)
    @test x(0.1,0.2) ≈ 0.1
    @test y(0.1,0.2) ≈ 0.2

    x,y=Fun(identity,d,20)
    @test x(0.1,0.2) ≈ 0.1
    @test y(0.1,0.2) ≈ 0.2


    # Boundary

    x,y=Fun(identity,∂(d),20)
    @test x(0.1,1.0) ≈ 0.1
    @test y(1.0,0.2) ≈ 0.2


    x,y=Fun(identity,∂(d))
    @test x(0.1,1.0) ≈ 0.1
    @test y(1.0,0.2) ≈ 0.2


    x,y=Fun(∂(d))
    @test x(0.1,1.0) ≈ 0.1
    @test y(1.0,0.2) ≈ 0.2
end

# test conversion between
dx=dy=Interval()
d=dx*dy

x,y=Fun(∂(d))
x,y=vec(x),vec(y)

g=[real(exp(x[1]-1im));0.0y[2];real(exp(x[3]+1im));real(exp(-1+1im*y[4]))]
B=[eye(dx)⊗ldirichlet(dy);ldirichlet(dx)⊗eye(dy);eye(dx)⊗rdirichlet(dy);rneumann(dx)⊗eye(dy)]

@test Fun(g[1],rangespace(B[1]))(-0.1,-1.0) ≈ g[1](-0.1,-1.0)
@test Fun(g[3],rangespace(B[3]))(-0.1,1.0)  ≈ g[3](-0.1,1.0)


A=ApproxFun.interlace([B;Δ])
g2=Fun([g;0.0],rangespace(A))

@test g2[1](-0.1,-1.0) ≈ g[1](-0.1,-1.0)
@test g2[3](-0.1,1.0)  ≈ g[3](-0.1,1.0)



S=WeightedJacobi(1,1)^2
L=Laplacian(S)
testbandedblockbandedoperator(L)



## Bug in Multiplication

dom = Interval(0.001, 1) * PeriodicInterval(-pi, pi)




r,r2 = Fun((r,t) -> [r;r^2], dom)

@test r(0.1,0.2) ≈ 0.1
@test r2(0.1,0.2) ≈ 0.1^2

sp = Space(dom)
Dr = Derivative(sp, [1,0])
Dθ = Derivative(sp, [0,1])
Mr = Multiplication(Fun( (r, θ) -> r, sp ), sp)
rDr = Mr * Dr

testbandedblockbandedoperator(rDr)



## Cheby * Interval
