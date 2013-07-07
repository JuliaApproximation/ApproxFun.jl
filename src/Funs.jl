module Funs
    using Base, Winston, NumericExtensions

export Fun,IFun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,plot,svfft




abstract AbstractFun
alternatingvector(n::Integer) = 2*mod([1:n],2)-1



include("Domains.jl")
include("IFun.jl")

##Fourier Routines
include("ShiftVector.jl")
include("FFun.jl")


##TODO: Add "FFun" for Fourier, following routine will decide
# whether input is IFun or FFun
Fun(x...)=apply(IFun,x)


end #module


