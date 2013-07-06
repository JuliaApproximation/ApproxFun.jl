using Funs

ef = IFun(exp);


@assert ef == -(-ef)
@assert ef == (ef-1) + 1


cf = IFun(cos); 

ecf = IFun(x->cos(x).*exp(x))

@assert abs(ef[.5]-exp(.5))<100eps()
@assert abs(ecf[.123456]-cos(.123456).*exp(.123456))<100eps()

r=2.*rand(100)-1;

@assert max(abs(ef[r]-exp(r)))<100eps()
@assert max(abs(ecf[r]-cos(r).*exp(r)))<100eps()


@assert norm((ecf-cf.*ef).coefficients)<100eps()

eocf = IFun(x->cos(x)./exp(x))

@assert max(abs((eocf-cf./ef).coefficients))<1000eps()


@assert norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


@assert norm((ef - diff(ef)).coefficients)<10E-11

@assert norm((ef - diff(cumsum(ef))).coefficients) < 10eps()
@assert norm((cf - diff(cumsum(cf))).coefficients) < 10eps()


ef2 = IFun(exp,[1,2]);
@assert norm((ef2 - diff(cumsum(ef2))).coefficients) < 10eps()


@assert abs(ef2[1.5]-exp(1.5))<100eps()

r=rand(100) + 1;

@assert max(abs(ef2[r]-exp(r)))<100eps()