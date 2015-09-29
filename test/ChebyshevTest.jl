using ApproxFun, Base.Test


ef = Fun(exp,Interval())


@test ef == -(-ef)
@test ef == (ef-1) + 1

ef = Fun(exp)


@test ef == -(-ef)
@test ef == (ef-1) + 1


cf = Fun(cos)

ecf = Fun(x->cos(x).*exp(x))
eocf = Fun(x->cos(x)./exp(x))

@test_approx_eq ef(.5) exp(.5)
@test_approx_eq ecf(.123456) cos(.123456).*exp(.123456)

r=2.*rand(100) .- 1

@test maximum(abs(ef(r)-exp(r)))<100eps()
@test maximum(abs(ecf(r)-cos(r).*exp(r)))<100eps()


@test norm((ecf-cf.*ef).coefficients)<100eps()



@test maximum(abs((eocf-cf./ef).coefficients))<1000eps()


@test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@test norm((ef - diff(ef)).coefficients)<10E-11

@test norm((ef - diff(cumsum(ef))).coefficients) < 20eps()
@test norm((cf - diff(cumsum(cf))).coefficients) < 20eps()

@test_approx_eq sum(ef)  2.3504023872876028

@test_approx_eq norm(ef)  1.90443178083307



##Check other domains


ef = Fun(exp,[1,2])
cf = Fun(cos,[1,2])

ecf = Fun(x->cos(x).*exp(x),[1,2])
eocf = Fun(x->cos(x)./exp(x),[1,2])


r=rand(100) .+ 1
x=1.5



@test_approx_eq ef(x) exp(x)



@test maximum(abs(ef(r)-exp(r)))<100eps()
@test maximum(abs(ecf(r)-cos(r).*exp(r)))<100eps()


@test norm((ecf-cf.*ef).coefficients)<100eps()


@test maximum(abs((eocf-cf./ef).coefficients))<1000eps()


@test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@test norm((ef - diff(ef)).coefficients)<10E-11

@test norm((ef - diff(cumsum(ef))).coefficients) < 10eps()
@test norm((cf - diff(cumsum(cf))).coefficients) < 10eps()

@test_approx_eq sum(ef) 4.670774270471604

@test_approx_eq norm(ef) 4.858451087240335


##Roots

f=Fun(x->sin(10(x-.1)))
@test norm(f(roots(f)))< 1000eps()


## ALiasing

f=Fun(x->cos(50acos(x)))
@test norm(f.coefficients-eye(length(f))[:,51])<100eps()


## Int values

@test_approx_eq Fun(x->2,10)(.1) 2
@test_approx_eq Fun(x->2)(.1) 2


@test_approx_eq Fun(Float64[],Chebyshev)([0.,1.]) [0.,0.]
@test_approx_eq Fun([],Chebyshev)(0.) 0.
@test_approx_eq Fun(x->[1.,0.])(0.) [1.,0.]
