using ApproxFun, Base.Test

@test norm(Fun(x->FFun(cos)[x])-Fun(cos)) <100eps()


@test norm(diff(FFun(cos))+FFun(sin)) < 100eps()

@test norm(diff(FFun(cos,Circle()))+FFun(sin,Circle()))<100eps()

f=Fun(exp,Circle());

@test norm(diff(f)-f)<100eps()
@test norm(integrate(f)+1-f)<100eps()