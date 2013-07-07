using Funs

@assert norm(IFun(x->FFun(cos,10)[x])-IFun(cos)) <10eps()


@assert norm(diff(f)+FFun(sin,10)) < 10eps()