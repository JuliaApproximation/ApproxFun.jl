@assert norm(IFun(x->FFun(cos)[x])-IFun(cos)) <100eps()


@assert norm(diff(FFun(cos))+FFun(sin)) < 100eps()

@assert norm(diff(FFun(cos,Circle()))+FFun(sin,Circle()))<100eps()

f=Fun(exp,Circle());

@assert norm(diff(f)-f)<100eps()
@assert norm(cumsum(f)+1-f)<100eps()