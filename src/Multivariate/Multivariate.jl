abstract MultivariateFun
abstract BivariateFun <:MultivariateFun


#implements coefficients/values/evaluate
Base.getindex(f::BivariateFun,x,y)=evaluate(f,x,y)
space(f::BivariateFun)=space(f,1)⊗space(f,2)
domain(f::BivariateFun)=domain(f,1)*domain(f,2)

differentiate(u::BivariateFun,i::Integer,j::Integer)=j==0?u:differentiate(differentiate(u,i),i,j-1)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)
lap(u::BivariateFun)=diff(u,1,2)+diff(u,2,2)
grad(u::BivariateFun)=[diff(u,1),diff(u,2)]

Base.chop(f::MultivariateFun)=chop(f,10eps())



include("VectorFun.jl")
include("TensorSpace.jl")
include("LowRankFun.jl")
include("TensorFun.jl")



Fun(f,S::TensorSpace,n...)=TensorFun(f,S,n...)
Fun(f,S::AbstractProductSpace,n...)=ProductFun(f,S,n...)
Fun(f,S::MultivariateDomain,n...)=Fun(f,Space(S),n...)
Fun(f,dx::Domain,dy::Domain)=Fun(f,dx*dy)
Fun(f,dx::Vector,dy::Vector)=Fun(f,Interval(dx),Interval(dx))

coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])
