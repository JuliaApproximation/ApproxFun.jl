using ApproxFun, Base.Test, DualNumbers


## Dual Numbers


f=Fun(exp,Segment(dual(1.0,1),dual(2.0)),20)
@test_approx_eq Fun(h->Fun(exp,Segment(1.0+h,2.0)).coefficients[1],0..1)'(0.) DualNumbers.epsilon(f.coefficients[1])

ud=let d=dual(0.0,1.0)..1.0
    B = ldirichlet(d)
    D = Derivative(d)
    a = Fun(exp,d)
    u = [B;D+dual(1.0,4.0)*a] \ [dual(1.0,2.0),0.0]
    u(0.5)
end

u0=let d=0.0..1.0
    B = ldirichlet(d)
    D = Derivative(d)
    a = Fun(exp,d)
    u = [B;D+a] \ [1.0,0.0]
    u(0.5)
end
h=0.00001
uh=let d=h..1.0
    B = ldirichlet(d)
    D = Derivative(d)
    a = Fun(exp,d)
    u = [B;D+(1+4h)*a] \ [1.0+2h,0.0]
    u(0.5)
end

@test_approx_eq_eps ud dual(u0,(uh-u0)/h) h

let d=0.0..1.0
    B = ldirichlet(d)
    D = Derivative(d)
    a = Fun(exp,d)
    u = [B;D+a] \ [dual(1.0,2.0),0.0]
    ur = [B;D+a] \ [1.0,0.0]
    ud = [B;D+a] \ [2.0,0.0]
    @test_approx_eq u(0.5)  dual(ur(0.5),ud(0.5))
end




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


## Sampling

ff=(x,y)->(x-y)^2*exp(-x^2/2-y^2/2)
ApproxFun.tensorizer(Chebyshev()^2)
f=Fun(ff,Domain(-4..4)^2)
r=ApproxFun.sample(f,5000)


#We can compare the histogram to the 1-point correlation
g=sum(f,1)/sum(f)
@test_approx_eq  g(0.1) 0.2004758624973169
