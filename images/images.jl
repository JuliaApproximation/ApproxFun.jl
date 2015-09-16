using PyCall
pygui(:tk)
using PyPlot

using ApproxFun
setplotter("PyPlot")
x = Fun(identity,[0.,10.])
f = sin(x^2)
g = cos(x)

h = f + g^2
r = roots(h)
rp = roots(differentiate(h))
ApproxFun.plot(h)               # using PyPlot
PyPlot.plot(r,h[r],"og",rp,h[rp],"or") # using PyPlot
xlabel("\$x\$");ylabel("\$h(x)\$");grid(true)
savefig("extrema.png",dpi=300)
clf()
println("First image done")


x = Fun(identity,[-1000.,200.])
d = domain(x)
D = Derivative(d)
B = dirichlet(d)
L = D^2 - x
u = [B,L] \ [airyai(d.a),airyai(d.b)]
ApproxFun.plot(u)						    # Requires Gadfly or PyPlot
xlabel("\$x\$");ylabel("\${\\rm Ai}(x)\$")
savefig("airy.png",dpi=300)
clf()
println("Second image done")


d = Interval([-π,π])
a = Fun(t-> 1+sin(cos(2t)),d)
D = Derivative(d)
L = D + a
f = Fun(t->exp(sin(10t)),d)
B = periodic(d,0)
uChebyshev = [B,L]\[0.,f]

d = PeriodicInterval([-π,π])
a = Fun(t-> 1+sin(cos(2t)),d)
D = Derivative(d)
L = D + a
f = Fun(t->exp(sin(10t)),d)
uFourier = L\f

length(uFourier)/length(uChebyshev),2/π
ApproxFun.plot(real(uFourier))						    # Requires Gadfly or PyPlot
xlim([d.a,d.b]);xlabel("\$t\$");ylabel("\$u(t)\$")
savefig("periodic.png",dpi=300)
clf()
println("Third image done")


f = abs(Fun(sin,[-5.,5.]))
d = domain(f)
x = ApproxFun.sample(f,10000)
ApproxFun.plot(f/sum(f))                           # Requires Gadfly or PyPlot
PyPlot.plt[:hist](x;normed=true,bins=[-5.:.1:5.],color="green")
xlim([first(d),last(d)]);ylim([0.0,0.18]);xlabel("\$x\$");ylabel("Density")
savefig("Sample.png",dpi=300)
clf()
println("Fourth image done")



## Nonlinear BVP

x=Fun()
u0=0.x

N=u->[u[-1.]-1.,u[1.]+0.5,0.001u''+6*(1-x^2)*u'+u^2-1.]
u=newton(N,u0)
ApproxFun.plot(u)
PyPlot.savefig("nbvp.png")
clf()
println("Fifth image done")
