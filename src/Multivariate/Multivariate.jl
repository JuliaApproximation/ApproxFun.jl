abstract MultivariateFun
#implements coefficients/values/evaluate
Base.getindex(f::MultivariateFun,x,y)=evaluate(f,x,y)



include("VectorFun.jl")
include("Fun2D.jl")
include("TensorSpace.jl")
include("TensorFun.jl")
