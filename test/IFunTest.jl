using Funs

ef = IFun(exp);


@assert ef == -(-ef)
@assert ef == (ef-1) + 1


cf = IFun(cos); 

ecf = IFun(x->cos(x).*exp(x))


@assert norm((ecf-cf.*ef).coefficients)<10eps()

eocf = IFun(x->cos(x)./exp(x))

@assert max(abs((eocf-cf./ef).coefficients))<1000eps()


@assert norm(((ef/3).*(3/ef)-1).coefficients)<1000eps()


@assert norm((ef - diff(ef)).coefficients)<10E-11