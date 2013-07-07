using Funs

@assert norm(IFun(x->FFun(cos,10)[x])-IFun(cos)) <10eps()


@assert norm(diff(FFun(cos,10))+FFun(sin,10)) < 10eps()

@assert norm(diff(FFun(cos,Circle(),10))+FFun(sin,Circle(),20))<10eps()

f=Fun(exp,Circle());

@assert norm(diff(f)-f)<10eps()
@assert norm(cumsum(f)+1-f)<10eps()