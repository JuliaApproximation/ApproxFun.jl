abstract MultivariateFun
abstract BivariateFun <:MultivariateFun


#implements coefficients/values/evaluate
Base.getindex(f::BivariateFun,x,y)=evaluate(f,x,y)
space(f::BivariateFun)=space(f,1)⊗space(f,2)
domain(f::BivariateFun)=domain(f,1)*domain(f,2)

differentiate(u::BivariateFun,i::Integer,j::Integer)=j==0?u:differentiate(differentiate(u,i),i,j-1)
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

arglength(f)=length(Base.uncompressed_ast(f.code.def).args[1])

function Fun(f::Function,d::BivariateSpace)
    if f==zero
        zeros(d)
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        Fun(LowRankFun(f,d))
    else
        Fun(LowRankFun((x,y)->f((x,y)),d))
    end
end


function Fun(f::Function,d::BivariateSpace,n::Integer)
    if (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        defaultFun(x->f(x...),d,n)
    else
        defaultFun(f,d,n)
    end
end

function Fun(f::Function)
    if (isgeneric(f)&&applicable(f,0)) || (!isgeneric(f)&&arglength(f)==1)
        # check for tuple
        try
            f(0)
        catch ex
            if isa(ex,BoundsError)
                # assume its a tuple
                return Fun(f,Interval()^2)
            else
                throw(ex)
            end
        end

        Fun(f,Interval())
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
            Fun(f,Interval()^2)
    else
        error("Function not defined on interval or square")
    end
end

Fun(f::ProductFun)=Fun(fromtensor(coefficients(f)),space(f))
Fun(f::ProductFun,sp::TensorSpace)=Fun(ProductFun(f,sp))
Fun(f::LowRankFun)=Fun(ProductFun(f))
Fun(f::LowRankFun,sp::TensorSpace)=Fun(ProductFun(f),sp)

coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])



points(f::BivariateFun,k...)=points(space(f),size(f,1),size(f,2),k...)


for OP in (:+,:-)
    @eval begin
        $OP(f::Fun,g::MultivariateFun)=$OP(ProductFun(f),g)
        $OP(f::MultivariateFun,g::Fun)=$OP(f,ProductFun(g))
    end
end
