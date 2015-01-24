abstract MultivariateFun
abstract BivariateFun <:MultivariateFun


#implements coefficients/values/evaluate
Base.getindex(f::BivariateFun,x,y)=evaluate(f,x,y)
space(f::BivariateFun)=space(f,1)⊗space(f,2)
domain(f::BivariateFun)=domain(f,1)*domain(f,2)

differentiate(u::BivariateFun,i::Integer,j::Integer)=j==0?u:differentiate(differentiate(u,i),i,j-1)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)
lap(u::BivariateFun)=differentiate(u,1,2)+differentiate(u,2,2)
grad(u::BivariateFun)=[differentiate(u,1),differentiate(u,2)]

Base.chop(f::MultivariateFun)=chop(f,10eps())



include("VectorFun.jl")
include("TensorSpace.jl")
include("LowRankFun.jl")
include("ProductFun.jl")


Fun(f,S::AbstractProductSpace,n...)=ProductFun(f,S,n...)
Fun(f,S::MultivariateDomain,n...)=Fun(f,Space(S),n...)
Fun(f,dx::Domain,dy::Domain)=Fun(f,dx*dy)
Fun(f,dx::Vector,dy::Vector)=Fun(f,Interval(dx),Interval(dx))


function Fun(f::Function)
    try
        Fun(f,Interval())
    catch ex #TODO only catch errors for wrong number of arguments
#    	warn("Got $(ex) when assuming 1-arity of $f")
#    	try
        	Fun(f,Interval(),Interval())
#         catch ex
#        	warn("Got $(ex) when assuming 2-arity of $f")
#         	error("Could not construct function")
#         end
    end
end


coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])


Base.zeros(sp::Union(MultivariateFunctionSpace,MultivariateDomain))=Fun(zeros(1,1),sp)
Base.zeros{T}(::Type{T},sp::Union(MultivariateFunctionSpace,MultivariateDomain))=Fun(zeros(T,1,1),sp)