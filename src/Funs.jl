module Funs
    using Base, Winston, NumericExtensions

export Fun,IFun,Interval,evaluate,values,points,chebyshev_transform
export pad!,pad,sample,chop!,complexroots,roots,plot




abstract AbstractFun
alternatingvector(n::Integer) = 2*mod([1:n],2)-1


include("Domains.jl");
include("IFun.jl");


##TODO: Add "FFun" for Fourier, following routine will decide
# whether input is IFun or FFun
Fun(x...)=apply(IFun,x)


end #module


