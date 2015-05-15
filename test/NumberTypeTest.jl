using ApproxFun, Base.Test

single_sin = Fun(sin,Interval(0.f0,1.f0))
double_sin = Fun(sin,Interval(0.,1.))
big_sin = Fun(sin,Interval(parse(BigFloat,"0.0"), parse(BigFloat,"1.0")))

@test length(single_sin) <= length(double_sin)
@test length(double_sin) <= length(big_sin)

@test eltype(coefficients(single_sin)) == Float32
@test eltype(coefficients(double_sin)) == Float64
@test eltype(coefficients(big_sin)) == BigFloat


single_double_err = coefficients(single_sin-double_sin)[1:length(single_sin)]
@test norm(single_double_err) < 10eps(Float32)

single_double_err = coefficients(double_sin-big_sin)[1:length(double_sin)]
@test norm(single_double_err) < 10eps(Float64)



