using ApproxFun, Base.Test

u0   = TensorFun((x,y)->cos(x)+sin(y) +exp(-50x.^2-40(y-.1).^2)+.5exp(-30(x+.5).^2-40(y+.2).^2))


@test values(u0)-values(u0|>Fun2D)|>norm < 1000eps()
@test chebyshevtransform(values(u0))-coefficients(u0)|>norm < 100eps()

##TODO: need to do adaptive to get better accuracy
@test sin(u0)[.1,.2]-sin(u0[.1,.2])|>abs < 10e-6
