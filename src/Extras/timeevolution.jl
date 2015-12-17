export timeevolution, BDF2,BDF4,BDF22

#
# These formulæ, appearing in Eq. (2.5) of:
#
# A.-K. Kassam and L. N. Trefethen, Fourth-order time-stepping for stiff PDEs, SIAM J. Sci. Comput., 26:1214--1233, 2005,
#
# are derived to implement ETDRK4 in double precision without numerical instability from cancellation.
#

expα(x) = (exp(x)*(4-3x+x^2)-4-x)/x^3
expβ(x) = (exp(x)*(x-2)+x+2)/x^3
expγ(x) = (exp(x)*(4-x)-4-3x-x^2)/x^3

expα_taylor(x::Float64) = Base.Math.@horner(x,1/6,1/6,3/40,1/45,5/1008,1/1120,7/51840,1/56700,1/492800,1/4790016,11/566092800,1/605404800,13/100590336000,1/106748928000,1/1580833013760,1/25009272288000,17/7155594141696000,1/7508956815360000)
expβ_taylor(x::Float64) = Base.Math.@horner(x,1/6,1/12,1/40,1/180,1/1008,1/6720,1/51840,1/453600,1/4435200,1/47900160,1/566092800,1/7264857600,1/100590336000,1/1494484992000,1/23712495206400,1/400148356608000,1/7155594141696000,1/135161222676480000)
expγ_taylor(x::Float64) = Base.Math.@horner(x,1/6,0/1,-1/120,-1/360,-1/1680,-1/10080,-1/72576,-1/604800,-1/5702400,-1/59875200,-1/691891200,-1/8717829120,-1/118879488000,-1/1743565824000,-1/27360571392000,-1/457312407552000,-1/8109673360588800)

expα(x::Float64) = abs2(x) < (17/16)^2 ? expα_taylor(x) : (exp(x)*(4-3x+x^2)-4-x)/x^3
expβ(x::Float64) = abs2(x) < (19/16)^2 ? expβ_taylor(x) : (exp(x)*(x-2)+x+2)/x^3
expγ(x::Float64) = abs2(x) < (15/16)^2 ? expγ_taylor(x) : (exp(x)*(4-x)-4-3x-x^2)/x^3

@vectorize_1arg Number expα
@vectorize_1arg Number expβ
@vectorize_1arg Number expγ


## Multivariate


function RK(L,y,h)
    k1=L(y)
    k2=L(y+.5h*k1)
    k3=L(y+.5h*k2)
    k4=L(y+h*k3)

    y+h*(k1+2k2+2k3+k4)/6.
end

function BDF2(B,A::BandedOperator,g::Function,bcs,u0,h,m,glp,tol=1000eps())
    SBDF2 = [B;I-2.0/3.0*h*A]

    u1=u0
    u2=chop(RK(g,u1,h),tol)
    u2,u1  = chop(SBDF2\[bcs;1/3.0*(4u2-u1)],tol),u2
    push!(glp,u2)
    yield()

    for k=1:m
        u2,u1 = chop(RK(g,u2,h),tol),u2
        u2,u1  = chop(SBDF2\[bcs;1/3.0*(4u2-u1)],tol),u2
        push!(glp,u2)
        yield()
    end

    u2
end




function BDF4(B::Vector,op::BandedOperator,bcs::Vector,uin::MultivariateFun,h::Real,m::Integer,glp)
    nt=size(uin,2)
    d=domain(uin)
    SBE   = discretize([B;I-h*op],d,nt)            # backward euler for first 2 time steps
    SBDF2 = discretize([B;I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF3 = discretize([B;I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF4 = discretize([B;I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps

    u1=uin
    u2=SBE\[bcs;u1]
    push!(glp,pad(u2,80,80))
    yield()
    u3=SBDF2\[bcs;1/3.0*(4u2-u1)]
    push!(glp,pad(u3,80,80))#updates window
    yield()
    u4=SBDF3\[bcs;1/11.0*(18u3-9u2+2u1)]
    push!(glp,pad(u4,80,80))#updates window
    yield()

    for k=1:m
        u4,u3,u2,u1  = SBDF4\[bcs;1/25.0*(48u4-36u3+16u2-3u1)],u4,u3,u2

        push!(glp,pad(u4,80,80))#updates window
        yield()
    end
    u4
end

BDF4(B::Vector,op::BandedOperator,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF4(B,op,zeros(length(B)),uin,h,m,glp)


function BDF22(B::Vector,op::BandedOperator,bcs::Vector,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,m::Integer,glp)
    nt=size(uin[1],2)
    SBE  = discretize([B;I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B;I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps

    u1,u2=uin
    u3 =SBE\[bcs;2u2-u1]
    push!(glp,pad(u3,80,80))
    yield()
    u4 =SBE\[bcs;2u3-u2]
    push!(glp,pad(u4,80,80))
    yield()

    for k=1:m
        u4,u3,u2,u1  = SBDF\[bcs;1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2

        push!(glp,pad(u4,80,80)) #updates window
        yield()
    end
    u4
end

function BDF22(B::Vector,op::BandedOperator,g::Function,bcs::Vector,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,m::Integer,glp)
    nt=size(uin[1],2)
    SBE  = discretize([B;I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B;I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps

    u1,u2=uin
    u3 =SBE\[bcs;2u2-u1]
    push!(glp,pad(u3,80,80))
    yield()
    u4,u3,u2,u1=u3,u2,u1,u1
    for k=1:m
        u4,u3,u2,u1 =2u4 - u3 +h^2*g(u4),u4,u3,u2
        u4,u3,u2,u1  = SBDF\[bcs;1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
        push!(glp,pad(u4,80,80))
        yield()
    end
    u4
end


BDF22(B::Vector,op::BandedOperator,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF22(B,op,zeros(length(B)),(uin,uin),h,m,glp)
BDF22(B::Vector,op::BandedOperator,g::Function,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF22(B,op,g,zeros(length(B)),(uin,uin),h,m,glp)


## GLPlot routines


#u_t = op*u
function timeevolution(B::Vector,op,bcs::Vector,uin::MultivariateFun,h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin,2)
    d=domain(uin)
    SBE   = discretize([B;I-h*op],d,nt)            # backward euler for first 2 time steps
    SBDF2 = discretize([B;I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF3 = discretize([B;I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF4 = discretize([B;I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps

    u1=chop(uin,1000eps())
    u2=chop(SBE\[bcs;u1],1000eps())
    plot(pad(u2,80,80),glp...)#updates window
    u3=chop(SBDF2\[bcs;1/3.0*(4u2-u1)],1000eps())
    plot(pad(u3,80,80),glp...)#updates window
    u4=chop(SBDF3\[bcs;1/11.0*(18u3-9u2+2u1)],1000eps())
    plot(pad(u4,80,80),glp...)#updates window

    for k=1:m
        u4,u3,u2,u1  = chop(SBDF4\[bcs;1/25.0*(48u4-36u3+16u2-3u1)],1000eps()),u4,u3,u2

        plot(pad(u4,80,80),glp...)#updates window
    end
    u4
end

function timeevolution(B::Vector,op,bcs::Vector,uin::MultivariateFun,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution(B,op,bcs,uin,h,m,plot(pad(uin,80,80)))
end

timeevolution(B::Vector,op,uin::MultivariateFun,h::Real,dat...)=timeevolution(B,op,zeros(length(B)),uin,h,dat...)
timeevolution(B::BandedOperator,dat...)=timeevolution([B],dat...)


timeevolution(B::Vector,op,bcs::Vector,uin::Fun,dat...)=timeevolution(B,op,bcs,ProductFun(uin),dat...)
timeevolution(B::Vector,op,uin::Fun,dat...)=timeevolution(B,op,ProductFun(uin),dat...)

timeevolution(B::Operator,dat...)=timeevolution([B],dat...)


#u_tt = op*u
function timeevolution2(B::Vector,op,bcs::Vector,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin[1],2)
    SBE  = discretize([B;I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B;I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps

    u1,u2=chop(uin[1],1000eps()),chop(uin[2],1000eps())
    u3 =chop(SBE\[bcs;2u2-u1],1000eps())
    u4 =chop(SBE\[bcs;2u3-u2],1000eps())

    for k=1:m
        u4,u3,u2,u1  = chop(SBDF\[bcs;1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2

        plot(pad(u4,80,80),glp...)#updates window
    end
    u4
end

#u_tt = op*u
function timeevolution2(B::Vector,op,g::Function,bcs::Vector,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin[1],2)
    SBE  = discretize([B;I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B;I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps

    u1,u2=chop(uin[1],1000eps()),chop(uin[2],1000eps())
    u3 =chop(SBE\[bcs;2u2-u1],1000eps())
    u4 =chop(2u3 - u2 +h^2*g(u3),1000eps())
    u4,u3,u2,u1  = chop(SBDF\[bcs;1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2
    plot(pad(u4,80,80),glp...)#updates window

    for k=1:m
        u4,u3,u2,u1  = chop(2u4 - u3 +h^2*g(u4),1000eps()),u4,u3,u2
        u4,u3,u2,u1  = chop(SBDF\[bcs;1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2
        plot(pad(u4,80,80),glp...)#updates window
    end
    u4
end

function timeevolution2(B::Vector,op,uin::Tuple{MultivariateFun,MultivariateFun},bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution2(B,op,bcs,uin,h,m,plot(pad(uin[end],80,80)))
end

function timeevolution2(B::Vector,op,g::Function,uin::Tuple{MultivariateFun,MultivariateFun},bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution2(B,op,g,bcs,uin,h,m,plot(pad(uin[end],80,80)))
end

timeevolution2(B::Vector,op,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,dat...)=timeevolution2(B,op,uin,zeros(length(B)),h,dat...)
timeevolution2(B::Vector,op,uin::MultivariateFun,dat...)=timeevolution2(B,op,(uin,uin),dat...)
timeevolution2(B::Operator,dat...)=timeevolution2([B],dat...)

timeevolution2(B::Vector,op,g::Function,uin::Tuple{MultivariateFun,MultivariateFun},h::Real,dat...)=timeevolution2(B,op,g,uin,zeros(length(B)),h,dat...)
timeevolution2(B::Vector,op,g::Function,uin::MultivariateFun,dat...)=timeevolution2(B,op,g,(uin,uin),dat...)


timeevolution2(B::Vector,op,g::Function,uin::Tuple{Fun,Fun},dat...)=timeevolution2(B,op,g,(ProductFun(uin[1]),ProductFun(uin[2])),dat...)
timeevolution2(B::Vector,op,g::Function,uin::Fun,dat...)=timeevolution2(B,op,g,ProductFun(uin),dat...)
timeevolution2(B::Vector,op,uin::Fun,dat...)=timeevolution2(B,op,ProductFun(uin),dat...)

timeevolution(o::Integer,dat...)=o==2?timeevolution2(dat...):timeevolution(dat...)
