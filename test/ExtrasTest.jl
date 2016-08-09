using ApproxFun, Base.Test


## Dual Numbers


f=Fun(exp,Interval(DualNumbers.dual(1.0,1),DualNumbers.dual(2.0)),20)
@test_approx_eq Fun(h->Fun(exp,Interval(1.0+h,2.0)).coefficients[1],[0.,.1])'(0.) DualNumbers.epsilon(f.coefficients[1])



## Eig test #336

d = Interval(0.,π)
A=Derivative(d)^2
λ=eigvals([dirichlet(d);A],100)
@test_approx_eq sort(λ)[end-5:end] -(-6:-1).^2


F = x->x.^8
d = Interval(0.0,1.0)
f = Fun(F,d)
ginf = Fun(x->exp(-x),d)
gp = ginf'
Af = Fun(x->x+f(x),d)
transport_ = Fun(x-> x - 1,d)
damping = Fun(x-> 1 - f(x),d)
A = transport_*Derivative(d) + damping
P = -DefiniteIntegral(Chebyshev(d))[LowRankFun((x,y)->gp(x)*(y+f(y)),d^2)];
λ,V = ApproxFun.eigs([A],100)
@test norm(sort(real(filter(x->isreal(x),λ)))[1:5]-(0:4)) ≤ 100000eps()

λ,V = ApproxFun.eigs([A+P],100)
@test_approx_eq_eps sort(real(filter(x->isreal(x),λ)))[5] 3.93759261234502 1E-3
