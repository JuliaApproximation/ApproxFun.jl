using ApproxFun, Base.Test


x=Fun(identity)
f=(exp(x)/(sqrt(1-x)*sqrt(x+1)))
@test_approx_eq (Hilbert(f|>space,0)*f)[.1] (-0.8545003781055088)
@test_approx_eq (Hilbert(0)*f)[.1] (-0.8545003781055088)
@test_approx_eq (Hilbert()*f)[.1] 1.1404096104609646386

x=Fun(identity,[-1,2])
f=(exp(x)/(sqrt(2-x)*sqrt(x+1)))
@test_approx_eq (Hilbert(f|>space,0)*f)[.1] 0.49127801561694168644
@test_approx_eq (Hilbert(0)*f)[.1] 0.49127801561694168644
@test_approx_eq (Hilbert()*f)[.1] 1.6649936695644078289



x=Fun(identity)
f=(exp(x)*(sqrt(1-x)*sqrt(x+1)))
@test_approx_eq (Hilbert()*f)[.1] 0.43723982258866913063

x=Fun(identity,[-1,2])
f=(exp(x)*(sqrt(2-x)*sqrt(x+1)))
@test_approx_eq (Hilbert()*f)[.1] 2.1380903070701673244