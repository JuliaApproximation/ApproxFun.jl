
ef = IFun(exp);


@assert ef == -(-ef)
@assert ef == (ef-1) + 1


cf = IFun(cos); 

ecf = IFun(x->cos(x).*exp(x))
eocf = IFun(x->cos(x)./exp(x))

@assert abs(ef[.5]-exp(.5))<100eps()
@assert abs(ecf[.123456]-cos(.123456).*exp(.123456))<100eps()

r=2.*rand(100)-1;

@assert max(abs(ef[r]-exp(r)))<100eps()
@assert max(abs(ecf[r]-cos(r).*exp(r)))<100eps()


@assert norm((ecf-cf.*ef).coefficients)<100eps()



@assert max(abs((eocf-cf./ef).coefficients))<1000eps()


@assert norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@assert norm((ef - diff(ef)).coefficients)<10E-11

@assert norm((ef - diff(cumsum(ef))).coefficients) < 10eps()
@assert norm((cf - diff(cumsum(cf))).coefficients) < 10eps()

@assert abs(sum(ef) - 2.3504023872876028) < 10eps()

@assert abs(norm(ef) - 1.90443178083307) < 10eps()



##Check other domains


ef = IFun(exp,[1,2]);
cf = IFun(cos,[1,2]); 

ecf = IFun(x->cos(x).*exp(x),[1,2])
eocf = IFun(x->cos(x)./exp(x),[1,2])


r=rand(100) + 1;
x=1.5;



@assert abs(ef[x]-exp(x))<100eps()



@assert max(abs(ef[r]-exp(r)))<100eps()
@assert max(abs(ecf[r]-cos(r).*exp(r)))<100eps()


@assert norm((ecf-cf.*ef).coefficients)<100eps()


@assert max(abs((eocf-cf./ef).coefficients))<1000eps()


@assert norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


## Diff and cumsum


@assert norm((ef - diff(ef)).coefficients)<10E-11

@assert norm((ef - diff(cumsum(ef))).coefficients) < 10eps()
@assert norm((cf - diff(cumsum(cf))).coefficients) < 10eps()

@assert abs(sum(ef) - 4.670774270471604) < 10eps()

@assert abs(norm(ef) - 4.858451087240335) < 10eps()

