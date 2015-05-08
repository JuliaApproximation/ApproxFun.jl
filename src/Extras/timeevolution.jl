export timeevolution, BDF2,BDF4,BDF22



## Multivariate 


function RK(L,y,h)
    k1=L(y)
    k2=L(y+.5h*k1)
    k3=L(y+.5h*k2)
    k4=L(y+h*k3)
    
    y+h*(k1+2k2+2k3+k4)/6.
end

function BDF2(B,A::BandedOperator,g::Function,bcs,u0,h,m,glp,tol=1000eps())
    SBDF2 = [B,I-2.0/3.0*h*A]

    u1=u0
    u2=chop(RK(g,u1,h),tol)
    u2,u1  = chop(SBDF2\[bcs,1/3.0*(4u2-u1)],tol),u2
    push!(glp,u2)

    for k=1:m
        u2,u1 = chop(RK(g,u2,h),tol),u2    
        u2,u1  = chop(SBDF2\[bcs,1/3.0*(4u2-u1)],tol),u2
        push!(glp,u2)
    end    

    u2
end




function BDF4(B::Vector,op::PDEOperator,bcs::Vector,uin::MultivariateFun,h::Real,m::Integer,glp)
    nt=size(uin,2)
    d=domain(uin)
    SBE   = discretize([B,I-h*op],d,nt)            # backward euler for first 2 time steps
    SBDF2 = discretize([B,I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF3 = discretize([B,I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps    
    SBDF4 = discretize([B,I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps        
    
    u1=uin
    u2=SBE\[bcs,u1]
    push!(glp,pad(u2,80,80))
    u3=SBDF2\[bcs,1/3.0*(4u2-u1)]
    push!(glp,pad(u3,80,80))#updates window        
    u4=SBDF3\[bcs,1/11.0*(18u3-9u2+2u1)]
    push!(glp,pad(u4,80,80))#updates window    
    
    for k=1:m
        u4,u3,u2,u1  = SBDF4\[bcs,1/25.0*(48u4-36u3+16u2-3u1)],u4,u3,u2
 
        push!(glp,pad(u4,80,80))#updates window
    end    
    u4
end

BDF4(B::Vector,op::PDEOperator,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF4(B,op,zeros(length(B)),uin,h,m,glp)


function BDF22(B::Vector,op::PDEOperator,bcs::Vector,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,m::Integer,glp)
    nt=size(uin[1],2)
    SBE  = discretize([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    push!(glp,pad(u3,80,80))
    u4 =SBE\[bcs,2u3-u2]
    push!(glp,pad(u4,80,80))
    
    for k=1:m
        u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
 
        push!(glp,pad(u4,80,80)) #updates window
    end    
    u4
end

function BDF22(B::Vector,op::PDEOperator,g::Function,bcs::Vector,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,m::Integer,glp)
    nt=size(uin[1],2)
    SBE  = discretize([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    push!(glp,pad(u3,80,80))
    u4,u3,u2,u1=u3,u2,u1,u1
    for k=1:m
        u4,u3,u2,u1 =2u4 - u3 +h^2*g(u4),u4,u3,u2
        u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2    
        push!(glp,pad(u4,80,80))
    end   
    u4 
end


BDF22(B::Vector,op::PDEOperator,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF22(B,op,zeros(length(B)),(uin,uin),h,m,glp)
BDF22(B::Vector,op::PDEOperator,g::Function,uin::MultivariateFun,h::Real,m::Integer,glp)=BDF22(B,op,g,zeros(length(B)),(uin,uin),h,m,glp)


## GLPlot routines


#u_t = op*u
function timeevolution(B::Vector,op,bcs::Vector,uin::MultivariateFun,h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin,2)
    d=domain(uin)
    SBE   = discretize([B,I-h*op],d,nt)            # backward euler for first 2 time steps
    SBDF2 = discretize([B,I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
    SBDF3 = discretize([B,I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps    
    SBDF4 = discretize([B,I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps        
    
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
function timeevolution2(B::Vector,op,bcs::Vector,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin[1],2)
    SBE  = discretize([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    u4 =SBE\[bcs,2u3-u2]

    for k=1:m
        u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
 
        plot(pad(u4,80,80),glp...)#updates window
    end    
end

#u_tt = op*u
function timeevolution2(B::Vector,op,g::Function,bcs::Vector,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    nt=size(uin[1],2)
    SBE  = discretize([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    SBDF = discretize([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    u4 =2u3 - u2 +h^2*g(u3)
    u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
    plot(pad(u4,80,80),glp...)#updates window                

    for k=1:m
        u4,u3,u2,u1  = chop(2u4 - u3 +h^2*g(u4),1000eps()),u4,u3,u2    
        u4,u3,u2,u1  = chop(SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2
        plot(pad(u4,80,80),glp...)#updates window               
    end    
end

function timeevolution2(B::Vector,op,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution2(B,op,bcs,uin,h,m,plot(pad(uin[end],80,80)))
end

function timeevolution2(B::Vector,op,g::Function,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),bcs::Vector,h::Real,m=5000)
    require("GLPlot")
    setplotter("GLPlot")
    timeevolution2(B,op,g,bcs,uin,h,m,plot(pad(uin[end],80,80)))
end

timeevolution2(B::Vector,op,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,dat...)=timeevolution2(B,op,uin,zeros(length(B)),h,dat...)
timeevolution2(B::Vector,op,uin::MultivariateFun,dat...)=timeevolution2(B,op,(uin,uin),dat...)
timeevolution2(B::PDEOperator,dat...)=timeevolution2([B],dat...)

timeevolution2(B::Vector,op,g::Function,uin::@compat(Tuple{MultivariateFun,MultivariateFun}),h::Real,dat...)=timeevolution2(B,op,g,uin,zeros(length(B)),h,dat...)
timeevolution2(B::Vector,op,g::Function,uin::MultivariateFun,dat...)=timeevolution2(B,op,g,(uin,uin),dat...)


timeevolution(o::Integer,dat...)=o==2?timeevolution2(dat...):timeevolution(dat...)

