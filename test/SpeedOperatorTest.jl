using ApproxFun, Base.Test


d=Interval(-1000.,5.)
x=Fun(identity,d)


A=Derivative(d)^2-x

@time A[1:12500,1:12505]  # 0.002191
@time A.ops[1][1:12500,1:12505]  # 0.000253
@time A.ops[2][1:12500,1:12505]  # 0.002507

# T=TimesOperator([A.ops[2].ops[1].op.ops...,A.ops[2].ops[2]])
# @time T[1:12500,1:12505]
# bandwidth(A.ops[2].ops[2],2)
@time a=A.ops[2].ops[1][1:12500,1:12504]  # 0.000770
@time b=A.ops[2].ops[2][1:12504,1:12505]  # 0.000117

@time A.ops[2].ops[1].op.ops[1][1:12500,1:12504]  # 0.000177
@time A.ops[2].ops[1].op.ops[2][1:12500,1:12504]  # 0.000167

@time A.ops[2].ops[2][1:12500,1:12504]  # 0.000091

#
# @profile for k=1:100
#     A.ops[2].ops[2][1:12500,1:12504]
# end
#
# Profile.print()
@time a*b  #0.0014
