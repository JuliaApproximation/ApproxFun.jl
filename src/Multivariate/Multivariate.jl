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


# Fun(f,S::MultivariateSpace,n...)=ProductFun(f,S,n...)
# Fun{T<:Number}(f::Number,S::MultivariateDomain{T})=Fun(f,Space(S))
# Fun{T<:Number}(f::Function,S::MultivariateDomain{T})=Fun(f,Space(S))
# Fun(f,S::MultivariateDomain,n...)=Fun(f,Space(S),n...)
# Fun{T<:Number}(f,dx::MultivariateDomain{T},dy::Domain)=Fun(f,dx*dy)
# Fun(f,dx::Domain,dy::Domain)=Fun(f,dx*dy)
# Fun(f,dx::Vector,dy::Vector)=Fun(f,Interval(dx),Interval(dx))


function Fun(f::Function)
    try
        f(0.)
        Fun(f,Interval())
    catch ex #TODO only catch errors for wrong number of arguments
#    	warn("Got $(ex) when assuming 1-arity of $f")
    	try
         	Fun(g,Interval()^2)
         catch ex
#        	warn("Got $(ex) when assuming 2-arity of $f")
         	error("Could not construct function")
         end
     end
end

Fun(f::ProductFun)=Fun(fromtensor(coefficients(f)),space(f))
function Fun(f::Function,S::BivariateSpace)
    try
        pt=checkpoints(S)[1]
        f(pt[1],pt[2])   # check if we can evaluate or need to dethread
        Fun(ProductFun(f,S))
    catch
        # assume it needs a tuple
        Fun(ProductFun((x,y)->f((x,y)),S))
    end 
end

coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])


Base.zeros(sp::Union(MultivariateSpace,MultivariateDomain))=Fun(zeros(1,1),sp)
Base.zeros{T<:Number}(::Type{T},sp::Union(MultivariateSpace,MultivariateDomain))=Fun(zeros(T,1,1),sp)



points(f::BivariateFun,k...)=points(space(f),size(f,1),size(f,2),k...)
