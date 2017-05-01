using ApproxFun, Base.Test


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
