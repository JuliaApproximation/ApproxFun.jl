abstract MultivariateFun

export Fun2D

#implements coefficients/values/evaluate
Base.getindex(f::MultivariateFun,x,y)=evaluate(f,x,y)
space(f::MultivariateFun)=space(f,1)⊗space(f,2)
domain(f::MultivariateFun)=domain(f,1)*domain(f,2)

differentiate(u::MultivariateFun,i::Integer,j::Integer)=j==0?u:differentiate(differentiate(u,i),i,j-1)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)
lap(u::MultivariateFun)=diff(u,1,2)+diff(u,2,2)
grad(u::MultivariateFun)=[diff(u,1),diff(u,2)]



include("VectorFun.jl")
include("TensorSpace.jl")
include("LowRankFun.jl")
include("TensorFun.jl")


Fun2D(f...)=LowRankFun(f...)