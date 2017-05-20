using ApproxFun, Base.Test



## broadcast

f=Fun(exp)
@test norm(exp.(f) - exp(f)) < 100eps()
@test norm(besselj.(1,f)-besselj(1,f)) < 1000eps()
@test atan2.(f,1)(0.1) ≈ atan2(f(0.1),1)
@test atan2.(f,f)(0.1) ≈ atan2(f(0.1),f(0.1))


x = Fun()
@test ≈(exp(x),exp.(x);atol=10eps())

f= Fun()
f .= exp.(x)
@test ≈(exp(x),f;atol=10eps())



f = Fun(Ultraspherical(1))
f .= exp.(x)
@test f(0.1) ≈ exp(0.1)


f = Fun()
@test_throws ArgumentError (f .= Fun(Line()))
