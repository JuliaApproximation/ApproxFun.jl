using PyCall
pygui(:tk)
using PyPlot
using Gadfly

using ApproxFun
setplotter("PyPlot")
x = Fun(identity,[0.,10.])
f = sin(x^2)
g = cos(x)

h = f + g^2
r = roots(h)
rp = roots(diff(h))
ApproxFun.plot(h)               # using PyPlot
PyPlot.plot(r,h[r],"og",rp,h[rp],"or") # using PyPlot
xlabel("\$x\$");ylabel("\$h(x)\$")
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
D = diff(d)
L = D + a
f = Fun(t->exp(sin(10t)),d)
B = periodic(d,0)
uChebyshev = [B,L]\[0.,f]

d = PeriodicInterval([-π,π])
a = FFun(t-> 1+sin(cos(2t)),d)
D = diff(d)
L = D + a
f = FFun(t->exp(sin(10t)),d)
uFourier = L\f

length(uFourier)/length(uChebyshev),2/π
ApproxFun.plot(real(uFourier))						    # Requires Gadfly or PyPlot
xlim([d.a,d.b]);xlabel("\$t\$");ylabel("\$u(t)\$")
savefig("periodic.png",dpi=300)
clf()
println("Third image done")


setplotter("Gadfly")
f = Fun(x->exp(-x^2),[-10,10])
x = ApproxFun.sample(f,10000)
ApproxFun.plot(f)             				# 2D plotting requires Gadfly or PyPlot
Gadfly.plot(x=x,Gadfly.Geom.histogram)
draw(PNG("Sample.png", 14inch, 7inch), Gadfly.plot(x=x,Gadfly.Geom.histogram))