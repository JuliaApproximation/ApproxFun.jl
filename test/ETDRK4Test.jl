using ApproxFun, Base.Test

import ApproxFun: expα, expβ, expγ

#
# This tests the numerical stability of the formulæ appearing in Eq. (2.5) of:
#
# A.-K. Kassam and L. N. Trefethen, Fourth-order time-stepping for stiff PDEs, SIAM J. Sci. Comput., 26:1214--1233, 2005.
#

x = logspace(-15,2,18)

@test norm(expα.(x)./expα.(big(x))-1,Inf) < eps()
@test norm(expβ.(x)./expβ.(big(x))-1,Inf) < eps()
@test norm(expγ.(x)./expγ.(big(x))-1,Inf) < eps()
