function timeevolutionsecondorder(B::Vector,op,bcs::Vector,uin::(MultivariateFun,MultivariateFun),h::Real,m::Integer,glp)
    require("GLPlot")
    setplotter("GLPlot")
    SBE  = schurfact([B,I-h^2*op],space(u0),100)            # backward euler for first 2 time steps
    SBDF = schurfact([B,I-4.0/9.0*h^2*op],space(u0),100)    # BDF formula for subsequent itme steps
    
    u1,u2=uin
    u3 =SBE\[bcs,2u2-u1]
    u4 =SBE\[bcs,2u3-u2]

    for k=1:m
        u4,u3,u2,u1  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],u4,u3,u2
 
        plot(pad(u4,50,50),glp...)#updates window
    end    
end

timeevolutionsecondorder(B::Vector,op,uin::(MultivariateFun,MultivariateFun),bcs::Vector,h::Real,m::Integer)=timeevolutionsecondorder(B,op,bcs,uin,h,m,plot(pad(uin[end],50,50)))

timeevolutionsecondorder(B::Vector,op,uin::(MultivariateFun,MultivariateFun),h::Real,dat...)=timeevolutionsecondorder(B,op,uin,zeros(length(B)),h,dat...)
timeevolutionsecondorder(B::Vector,op,uin::MultivariateFun,dat...)=timeevolutionsecondorder(B,op,(uin,uin),dat...)
timeevolutionsecondorder(B::PDEOperator,dat...)=timeevolutionsecondorder([B],dat...)