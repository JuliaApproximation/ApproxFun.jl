# ## System of nonlinear differential equations
# One can also solve a system of nonlinear ODEs with potentially nonlinear boundary conditions:

using ApproxFun
using LinearAlgebra

N(u1, u2) = [u1'(0) - 0.5*u1(0)*u2(0);
                u2'(0) + 1;
                u1(1) - 1;
                u2(1) - 1;
                u1'' + u1*u2;
                u2'' - u1*u2]

function nbvpsolver2()
    x = Fun(0..1)
    u10 = one(x)
    u20 = one(x)
    newton(N, [u10,u20])
end
u1,u2 = nbvpsolver2();
