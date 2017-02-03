using ApproxFun, Base.Test
    import ApproxFun: testspace

for d in (Interval(),Interval(1.,2.),Segment(1.0+im,2.0+2im))
    testspace(Chebyshev(d))
end


ef = Fun(exp,Interval())

@test ef == -(-ef)
@test ef == (ef-1) + 1

ef = Fun(exp)


@test ef == -(-ef)
@test ef == (ef-1) + 1


cf = Fun(cos)

ecf = Fun(x->cos(x).*exp(x))
eocf = Fun(x->cos(x)./exp(x))

@test ef(.5) ≈ exp(.5)
@test ecf(.123456) ≈ cos(.123456).*exp(.123456)

r=2.*rand(100) .- 1

@test maximum(abs(ef(r)-exp(r)))<100eps()
@test maximum(abs(ecf(r)-cos(r).*exp(r)))<100eps()


@test norm((ecf-cf.*ef).coefficients)<100eps()



@test maximum(abs((eocf-cf./ef).coefficients))<1000eps()


@test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@test norm((ef - ef').coefficients)<10E-11

@test norm((ef - cumsum(ef)').coefficients) < 20eps()
@test norm((cf - cumsum(cf)').coefficients) < 20eps()

@test sum(ef)  ≈ 2.3504023872876028

@test norm(ef)  ≈ 1.90443178083307



##Check other domains


ef = Fun(exp,1..2)
cf = Fun(cos,1..2)

ecf = Fun(x->cos(x).*exp(x),1..2)
eocf = Fun(x->cos(x)./exp(x),1..2)


r=rand(100) .+ 1
x=1.5



@test ef(x) ≈ exp(x)



@test maximum(abs(ef(r)-exp(r)))<100eps()
@test maximum(abs(ecf(r)-cos(r).*exp(r)))<100eps()


@test norm((ecf-cf.*ef).coefficients)<100eps()


@test maximum(abs((eocf-cf./ef).coefficients))<1000eps()


@test norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@test norm((ef - ef').coefficients)<10E-11

@test norm((ef - cumsum(ef)').coefficients) < 10eps()
@test norm((cf - cumsum(cf)').coefficients) < 10eps()

@test sum(ef) ≈ 4.670774270471604

@test norm(ef) ≈ 4.858451087240335


##Roots

f=Fun(x->sin(10(x-.1)))
@test norm(f(roots(f)))< 1000eps()


## ALiasing

f=Fun(x->cos(50acos(x)))
@test norm(f.coefficients-eye(ncoefficients(f))[:,51])<100eps()


## Int values

@test Fun(x->2,10)(.1) ≈ 2
@test Fun(x->2)(.1) ≈ 2


@test Fun(Chebyshev,Float64[])([0.,1.]) ≈ [0.,0.]
@test Fun(Chebyshev,[])(0.) ≈ 0.
@test Fun(x->[1.,0.])(0.) ≈ [1.,0.]



## broadcast

f=Fun(exp)
@test norm(exp.(f) - exp(f)) < 100eps()
@test norm(besselj.(1,f)-besselj(1,f)) < 100eps()
@test atan2.(f,1)(0.1) ≈ atan2(f(0.1),1)
@test atan2.(f,f)(0.1) ≈ atan2(f(0.1),f(0.1))




## Fixes #121

x = Fun(identity,0..10)
f = sin(x^2)
g = cos(x)
@test f(.1) ≈ sin(.1^2)

x = Fun(identity,0..100)
f = sin(x^2)
@test_approx_eq_eps f(.1) sin(.1^2) 1E-12


## Reverse


f=Fun(exp)
@test Fun(f,Chebyshev(1..(-1)))(0.1) ≈ f(0.1)
