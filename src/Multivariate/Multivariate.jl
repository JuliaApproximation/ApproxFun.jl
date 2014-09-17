abstract MultivariateFun

#implements coefficients/values/evaluate
Base.getindex(f::MultivariateFun,x,y)=evaluate(f,x,y)
space(f::MultivariateFun)=space(f,1)⊗space(f,2)
domain(f::MultivariateFun)=domain(f,1)*domain(f,2)

Base.diff(u::MultivariateFun,i::Integer,j::Integer)=j==0?u:diff(diff(u,i),i,j-1)
lap(u::MultivariateFun)=diff(u,1,2)+diff(u,2,2)
grad(u::MultivariateFun)=[diff(u,1),diff(u,2)]



include("VectorFun.jl")
include("TensorSpace.jl")
include("Fun2D.jl")
include("TensorFun.jl")
