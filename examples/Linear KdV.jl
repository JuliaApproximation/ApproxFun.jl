using ApproxFun
import ApproxFun: factor, component

L = 1.0
T = 0.1
d = Interval(0.,T) * Interval(0.,L)
Dx = Derivative(d,[0,1])
Dt = Derivative(d,[1,0])

α = 1.0

B = [ldirichlet() ⊗ I;
     I ⊗ ldirichlet();
     I ⊗ rdirichlet();
     (I ⊗ rneumann()) - α*(I ⊗ lneumann());]

q0 = Fun(x->x^2*(L-x)^2 , factor(d,2))
f0 = zeros(component(∂(d),1))
g0 = zeros(component(∂(d),3))

u = \([B;Dt + Dx^3] , [q0;f0;g0;0;0] ; tolerance=1E-5)
