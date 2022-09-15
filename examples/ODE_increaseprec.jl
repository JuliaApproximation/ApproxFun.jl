using ApproxFun
u = setprecision(1000) do
    d = BigFloat(0)..BigFloat(1)
    D = Derivative(d)
    [ldirichlet(); D-1] \ [1; 0]
end;
using Test #src
@test u(1) ≈ exp(BigFloat(1)) #src
