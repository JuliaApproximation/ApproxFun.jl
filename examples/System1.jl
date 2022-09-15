using ApproxFun
using LinearAlgebra

x = Fun();
B = Evaluation(Chebyshev(),-1);
A = [B      0;
     B*𝒟    0;
     0      B;
     𝒟^2-I  2I;
     I      𝒟+I];
u,v = A \ [0;0;0;exp(x);cos(x)];
