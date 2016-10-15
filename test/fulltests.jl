## This includes extra tests that are too time consuming for Travis


include("runtests.jl")


println("    Bessel tests")

@time for ν in (1.,0.5,2.,3.5)
    println("        ν = $ν")
    S=JacobiWeight(-ν,0.,Chebyshev([0.,1.]))
    D=Derivative(S)
    x=Fun(identity,domain(S))
    L=(x^2)*D^2+x*D+(x^2-ν^2);
    u=linsolve([rdirichlet(S);rneumann(S);L],[bessely(ν,1.),.5*(bessely(ν-1.,1.)-bessely(ν+1.,1.))];
                tolerance=1E-13)
    @test_approx_eq_eps u(.1) bessely(ν,.1) eps(100000.)*max(abs(u(.1)),1)
    u=linsolve([rdirichlet(S),rneumann(S),L],[besselj(ν,1.),.5*(besselj(ν-1.,1.)-besselj(ν+1.,1.))];
                tolerance=1E-13)
    @test_approx_eq_eps u(.1) besselj(ν,.1) eps(100000.)*max(abs(u(.1)),1)

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
