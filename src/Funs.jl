module Funs
    using Base, Winston

export IFun,Interval,evaluate,values,points,chebyshev_transform
export pad!,pad,sample




include("Domains.jl");
include("IFun.jl");


end #module


