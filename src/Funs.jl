module Funs
    using Base, Winston

export Fun,IFun,Interval,evaluate,values,points,chebyshev_transform
export pad!,pad,sample,chop!,complexroots,roots


abstract AbstractFun

include("Domains.jl");
include("IFun.jl");


end #module


