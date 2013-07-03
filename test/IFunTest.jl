

ef = IFun(exp);


@assert ef == -(-ef)
@assert ef == (ef-1) + 1


cf = IFun(cos); 

ecf = IFun(x->cos(x).*exp(x))


@assert norm((ecf-cf.*ef).coefficients)<10eps()

eocf = IFun(x->cos(x)./exp(x))

@assert norm((eocf-cf./ef).coefficients)<10eps()


@assert norm(((cf/3).*(3/cf)-1).coefficients)<100eps()