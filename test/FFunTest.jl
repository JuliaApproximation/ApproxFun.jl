using Funs

@assert norm(IFun(x->FFun(cos,10)[x])-IFun(cos)) <10eps()