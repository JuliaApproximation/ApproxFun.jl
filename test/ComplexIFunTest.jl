using ApproxFun, Base.Test

## Diff


f=Fun(x->exp(im.*x))

@test norm(f'-im*f) < 1000eps()


@test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()


@test norm(real(f)-Fun(cos)) < eps()

@test norm(real(f-Fun(cos))) < eps()


##Check other real domains


f=Fun(x->exp(im.*x),[1,2])


@test norm(f'-im*f) < 1000eps()


@test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()


@test norm(real(f)-Fun(cos,domain(f))) < eps()

@test norm(real(f-Fun(cos,domain(f)))) < eps()


##Check complex domains


f=Fun(x->exp(im.*x),[1im,2+.5im])

#@assert f([f.domain.a,f.domain.b])      ##TODO: Currently crashes


@test_approx_eq sum(f)  (0.5515167681675808 + 0.6202852564797062im)

@test_approx_eq f(1im) exp(im.*im)

@test_approx_eq f(1+.75im) exp(im.*(1+.75im))


@test norm(f'-im*f) < 1000eps()


@test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()
