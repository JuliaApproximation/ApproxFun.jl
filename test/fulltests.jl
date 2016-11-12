using ApproxFun, Base.Test
    import ApproxFun: testbandedoperator, testraggedbelowoperator, InterlaceOperator
## This includes extra tests that are too time consuming for Travis


include("runtests.jl")

println("    Full Operator tests")

@time for M in (Multiplication(Fun([1.],CosSpace()),CosSpace()),
                Multiplication(Fun([1.],CosSpace()),SinSpace()),
                Multiplication(Fun([1.],SinSpace()),SinSpace()),
                Multiplication(Fun([1.],SinSpace()),CosSpace()),
                Derivative(SinSpace()),Derivative(CosSpace()))
      testbandedoperator(M)
end

S=Chebyshev()
@time for io in (
        InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev());lneumann(S)]),
        InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos);lneumann(S)]),
        InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())]),
        InterlaceOperator([InterlaceOperator(dirichlet(S));Derivative(Chebyshev())+Fun(cos)]),
        InterlaceOperator([Derivative(Chebyshev());InterlaceOperator(dirichlet(S))]),
        InterlaceOperator([Derivative(Chebyshev())+Fun(cos);InterlaceOperator(dirichlet(S))]))
    testraggedbelowoperator(io)
end

## Newton iteration bug
S=Chebyshev([0.,7.])

ω=2π


N = u->Any[Fun(u(0.)-0.1);Fun(u(ω)-u(0.));Fun(u'(ω)-u'(0.));u''+u+u^3]

u=0.1Fun(cos,S)

D=Derivative(S)

Z=ApproxFun.ZeroOperator(ApproxFun.ConstantSpace())

A=ApproxFun.interlace([Z                      Evaluation(S,0);
                     u'(ω)    Evaluation(S,ω)-Evaluation(S,0);
                     u''(ω)   Evaluation(S,ω,1)-Evaluation(S,0,1);
                      0         D^2+I+3u^2])

@time testraggedbelowoperator(A)

# This is surprisingly slow...

DS=WeightedJacobi(0.1+1,0.2+1)
D=Derivative(DS)[2:end,:]

@time ApproxFun.testbandedoperator(D)


println("    Bessel tests")

@time for ν in (1.,0.5,2.,3.5)
    println("        ν = $ν")
    S=JacobiWeight(-ν,0.,Chebyshev([0.,1.]))
    D=Derivative(S)
    x=Fun(identity,domain(S))
    L=(x^2)*D^2+x*D+(x^2-ν^2);
    u=linsolve([rdirichlet(S);rneumann(S);L],[bessely(ν,1.),.5*(bessely(ν-1.,1.)-bessely(ν+1.,1.))];
                tolerance=1E-10)
    @test_approx_eq_eps u(.1) bessely(ν,.1) eps(1000000.)*max(abs(u(.1)),1)
    u=linsolve([rdirichlet(S),rneumann(S),L],[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.))];
                tolerance=1E-10)
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(1000000.)*max(abs(u(.1)),1)

    u=Fun(x->bessely(ν,x),S)
    @test_approx_eq_eps u(.1) bessely(ν,.1) eps(10000.)*max(abs(u(.1)),1)
    u=Fun(x->besselj(ν,x),S)
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(10000.)*max(abs(u(.1)),1)
end

@time for ν in (1.,0.5,0.123,3.5)
    println("        ν = $ν")
    S=JacobiWeight(ν,0.,Chebyshev([0.,1.]))
    D=Derivative(S)
    x=Fun(identity,domain(S))
    L=(x^2)*D^2+x*D+(x^2-ν^2);

    u=linsolve([rdirichlet(S),rneumann(S),L],[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.))];
                tolerance=1E-10)
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(1000000.)*max(abs(u(.1)),1)

    u=Fun(x->besselj(ν,x),S)
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(10000.)*max(abs(u(.1)),1)
end

println("Full Jacobi tests")

## Legendre conversions
testspace(Ultraspherical(1);haslineintegral=false)
testspace(Ultraspherical(2);haslineintegral=false)
@time testspace(Ultraspherical(1//2);haslineintegral=false,minpoints=2)  # minpoints is a tempory fix a bug

@test norm(Fun(exp,Ultraspherical(1//2))-Fun(exp,Jacobi(0,0))) < 100eps()

C=Conversion(Jacobi(0,0),Chebyshev())
@time testbandedbelowoperator(C)
@test norm(C*Fun(exp,Jacobi(0,0))  - Fun(exp)) < 100eps()


C=Conversion(Ultraspherical(1//2),Chebyshev())
@time testbandedbelowoperator(C)
@test norm(C*Fun(exp,Ultraspherical(1//2))  - Fun(exp)) < 100eps()



C=Conversion(Chebyshev(),Ultraspherical(1//2))
@time testbandedbelowoperator(C)
@test norm(C*Fun(exp)-Fun(exp,Legendre())) < 100eps()


C=Conversion(Chebyshev(),Jacobi(0,0))
@time testbandedbelowoperator(C)
@test norm(C*Fun(exp)  - Fun(exp,Jacobi(0,0))) < 100eps()


C=Conversion(Chebyshev(),Jacobi(1,1))
@time testbandedbelowoperator(C)
@test norm(C*Fun(exp) - Fun(exp,Jacobi(1,1))) < 100eps()


C=Conversion(Ultraspherical(1//2),Ultraspherical(1))
@time testbandedbelowoperator(C)

λ1 = ApproxFun.order(domainspace(C))
λ2 = ApproxFun.order(rangespace(C))

# test against version that doesn't use lgamma
Cex = Float64[(if j ≥ k && iseven(k-j)
        gamma(λ2)*(k-1+λ2)/(gamma(λ1)*gamma(λ1-λ2))*
            (gamma((j-k)/2+λ1-λ2)/gamma((j-k)/2+1))*
            (gamma((k+j-2)/2+λ1)/gamma((k+j-2)/2+λ2+1))
    else
        0.0
    end) for k=1:20,j=1:20]

@test norm(Cex - C[1:20,1:20]) < 100eps()

@test norm(C*Fun(exp,Ultraspherical(1//2))-Fun(exp,Ultraspherical(1))) < 100eps()

C=Conversion(Jacobi(0,0),Ultraspherical(1))
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Jacobi(0,0))-Fun(exp,Ultraspherical(1))) < 100eps()


C=Conversion(Ultraspherical(1),Jacobi(0,0))
testbandedbelowoperator(C)
@test norm(C*Fun(exp,Ultraspherical(1))-Fun(exp,Jacobi(0,0))) < 100eps()


println("Full multivariate tests")


## ProductFun
u0   = ProductFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>LowRankFun)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)(.1,.2)-sin(u0(.1,.2))|>abs < 10e-4


F = LowRankFun((x,y)->hankelh1(0,10abs(y-x)),Chebyshev([1.0,2.0]),Chebyshev([1.0im,2.0im]))

@test_approx_eq F(1.5,1.5im) hankelh1(0,10abs(1.5im-1.5))


## Periodic
f=LowRankFun((x,y)->cos(x)*sin(y),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f(.1,.2) cos(.1)*sin(.2)

f=LowRankFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval(),PeriodicInterval())
@test_approx_eq f(.1,.2) cos(cos(.1)+sin(.2))
@test norm(Float64[cos(cos(x)+sin(y)) for x=ApproxFun.vecpoints(f,1),y=ApproxFun.vecpoints(f,2)]-values(f))<10000eps()

f=ProductFun((x,y)->cos(cos(x)+sin(y)),PeriodicInterval()^2)
@test_approx_eq f(.1,.2) cos(cos(.1)+sin(.2))
x,y=points(f)
@test norm(Float64[cos(cos(x[k,j])+sin(y[k,j])) for k=1:size(f,1),j=1:size(f,2)]-values(f))<10000eps()

d=PeriodicInterval()^2
f=ProductFun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
@test (f.'-f|>coefficients|>norm)< 1000eps()

## Functional*Fun

d=Interval()
B=ldirichlet(d)
f=ProductFun((x,y)->cos(cos(x)*sin(y)),d^2)
@test norm(B*f-Fun(y->cos(cos(-1)*sin(y)),d))<20000eps()
@test norm(f*B-Fun(x->cos(cos(x)*sin(-1)),d))<20000eps()

## matrix

f=Fun((x,y)->[exp(x*cos(y));cos(x*sin(y));2],Interval()^2)
@test_approx_eq f(0.1,0.2) [exp(0.1*cos(0.2));cos(0.1*sin(0.2));2]

f=Fun((x,y)->[exp(x*cos(y)) cos(x*sin(y)); 2 1],Interval()^2)
@test_approx_eq f(0.1,0.2) [exp(0.1*cos(0.2)) cos(0.1*sin(0.2));2 1]


## Cauchy fun

f = Fun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),Line()^2)
@test_approx_eq f(0.1,0.2) 1/(2π*(0.1^2+0.2^2+1)^(3/2))


#TODO: improve tolerance
f = LowRankFun((x,y)->1/(2π*(x^2+y^2+1)^(3/2)),JacobiWeight(2.,2.,Line())^2)
@test_approx_eq_eps f(0.1,0.2) 1/(2π*(0.1^2+0.2^2+1)^(3/2)) 1E-4





println("Full PDE tests")

include("FullPDETest.jl")


println("Speed tests")
include("SpeedTest.jl")
include("SpeedODETest.jl")
include("SpeedPDETest.jl")



println("Example tests")
if isdir(Pkg.dir("GR")) || isdir(Pkg.dir("Plotly")) || isdir(Pkg.dir("PlotlyJS")) || isdir(Pkg.dir("PyPlot"))
    include("ExamplesTest.jl")
else
    warn("Unable to do Examples since Gadfly.jl is not installed")
end


println("Readme tests")
include("ReadmeTest.jl")
