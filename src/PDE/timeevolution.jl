export timeevolution

#u_t = op*u
function timeevolution(B::Vector,op,bcs::Vector,uin::MultivariateFun,h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin,2)
    d=domain(uin)
    SBE   = schurfact([B,I-h*op],d,nt)            # backward euler for first 2 time steps
    SBDF2 = schurfact([B,I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF3 = schurfact([B,I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps    
    SBDF4 = schurfact([B,I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps        
    
    u1=uin
    u2=SBE\[bcs,u1]
    plot(pad(u2,80,80),glp...)#updates window    
    u3=SBDF2\[bcs,1/3.0*(4u2-u1)]
    plot(pad(u3,80,80),glp...)#updates window        
    u4=SBDF3\[bcs,1/11.0*(18u3-9u2+2u1)]
    plot(pad(u4,80,80),glp...)#updates window    
    
    for k=1:m
        u4,u3,u2,u1  = SBDF4\[bcs,1/25.0*(48u4-36u3+16u2-3u1)],u4,u3,u2
 
        plot(pad(u4,80,80),glp...)#updates window
    end    
end

function timeevolution(B::Vector,op,uin::MultivariateFun,bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution(B,op,bcs,uin,h,m,plot(pad(uin,80,80)))
end

timeevolution(B::Vector,op,uin::MultivariateFun,h::Real,dat...)=timeevolution(B,op,uin,zeros(length(B)),h,dat...)
timeevolution(B::PDEOperator,dat...)=timeevolution([B],dat...)



#u_tt = op*u
function timeevolution2(B::Vector,op,bcs::Vector,uin::(MultivariateFun,MultivariateFun),h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin[1],2)
    SBE  = schurfact([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = schurfact([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    u4 =SBE\[bcs,2u3-u2]

    for k=1:m
        u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
 
        plot(pad(u4,80,80),glp...)#updates window
    end    
end

function timeevolution2(B::Vector,op,uin::(MultivariateFun,MultivariateFun),bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution2(B,op,bcs,uin,h,m,plot(pad(uin[end],80,80)))
end

timeevolution2(B::Vector,op,uin::(MultivariateFun,MultivariateFun),h::Real,dat...)=timeevolution2(B,op,uin,zeros(length(B)),h,dat...)
timeevolution2(B::Vector,op,uin::MultivariateFun,dat...)=timeevolution2(B,op,(uin,uin),dat...)
timeevolution2(B::PDEOperator,dat...)=timeevolution2([B],dat...)


timeevolution(o::Integer,dat...)=o==2?timeevolution2(dat...):timeevolution(dat...)