

ef = IFun(exp);


@assert ef == -(-ef)
@assert ef == (ef-1) + 1


cf = IFun(cos); 

ecf = IFun(x->cos(x).*exp(x))