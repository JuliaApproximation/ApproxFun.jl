using Plots,ApproxFun;pyplot()

x = Fun(identity,[0.,10.])
f = sin(x^2)
g = cos(x)

h = f + g^2
r = roots(h)
rp = roots(differentiate(h))
plot(h;legend=false,grid=false)
scatter!(r,h(r))
scatter!(rp,h(rp))
xlabel!("\$x\$")
ylabel!("\$h(x)\$")
xlims!(0,10)
PyPlot.savefig("extrema.png",dpi=300)
println("First image done")


x = Fun(identity,[-1000.,200.])
d = domain(x)
D = Derivative(d)
B = dirichlet(d)
L = D^2 - x
u = [B;L] \ [airyai(d.a),airyai(d.b)]
plot(u;legend=false,grid=false)						    # Requires Gadfly or PyPlot
xlabel!("\$x\$");ylabel!("\${\\rm Ai}(x)\$")
PyPlot.savefig("airy.png",dpi=300)

println("Second image done")


d = Interval([-π,π])
a = Fun(t-> 1+sin(cos(2t)),d)
D = Derivative(d)
L = D + a
f = Fun(t->exp(sin(10t)),d)
B = periodic(d,0)
uChebyshev = [B;L]\[0.,f]

d = PeriodicInterval([-π,π])
a = Fun(t-> 1+sin(cos(2t)),d)
D = Derivative(d)
L = D + a
f = Fun(t->exp(sin(10t)),d)
uFourier = L\f

length(uFourier)/length(uChebyshev),2/π
plot(real(uFourier);legend=false,grid=false)						    # Requires Gadfly or PyPlot
xlims!(d.a,d.b);xlabel!("\$t\$");ylabel!("\$u(t)\$")
PyPlot.savefig("periodic.png",dpi=300)

println("Third image done")


f = abs(Fun(sin,[-5.,5.]))
d = domain(f)
x = ApproxFun.sample(f,10000)
plot(f/sum(f);legend=false,grid=false)                           # Requires Gadfly or PyPlot
histogram!(x;normed=true,nbins=100)
xlims!(first(d),last(d));ylims!(0.0,0.18);xlabel!("\$x\$");ylabel!("Density")
PyPlot.savefig("Sample.png",dpi=300)

println("Fourth image done")



## Nonlinear BVP

x=Fun()
u0=0.x

N=u->[u(-1.)-1.,u(1.)+0.5,0.001u''+6*(1-x^2)*u'+u^2-1.]
u=newton(N,u0)
plot(u;legend=false,grid=false)
PyPlot.savefig("nbvp.png",dpi=300)

println("Fifth image done")


## Multivariate


d = Interval()^2                            # Defines a rectangle

u = [dirichlet(d);lap(d)+100I]\ones(4)      # First four entries of rhs are
                                            # boundary conditions
surface(u,grid=false)                                     # contour plot

PyPlot.savefig("helmholtz.png",dpi=300)
