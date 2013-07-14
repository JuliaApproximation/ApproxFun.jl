using ApproxFun

## Diff


f=Fun(x->exp(im.*x))

@assert norm(diff(f)-im*f) < 1000eps()


@assert norm(cumsum(f)+im*f-f.coefficients[1]*im) < 100eps()


@assert norm(real(f)-Fun(cos)) < eps()

@assert norm(real(f-Fun(cos))) < eps()


##Check other real domains


f=Fun(x->exp(im.*x),[1,2])


@assert norm(diff(f)-im*f) < 1000eps()


@assert norm(cumsum(f)+im*f-f.coefficients[1]*im) < 100eps()


@assert norm(real(f)-Fun(cos,f.domain)) < eps()

@assert norm(real(f-Fun(cos,f.domain))) < eps()


##Check complex domains


f=Fun(x->exp(im.*x),[1im,2+.5im])

@assert f[[f.domain.a,f.domain.b]]      ##TODO: Currently crashes


@assert abs(sum(f) - (0.5515167681675808 + 0.6202852564797062im)) < 10eps()

@assert abs(f[1im]-exp(im.*im)) <eps()

@assert abs(f[1+.75im]-exp(im.*(1+.75im))) < 100eps()


@assert norm(diff(f)-im*f) < 1000eps()


@assert norm(cumsum(f)+im*f-f.coefficients[1]*im) < 100eps()

