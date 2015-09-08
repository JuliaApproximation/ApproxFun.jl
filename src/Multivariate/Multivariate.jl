abstract MultivariateFun{T}
abstract BivariateFun{T} <: MultivariateFun{T}

export grad,lap,curl

#implements coefficients/values/evaluate
Base.getindex(f::BivariateFun,x,y)=evaluate(f,x,y)
space(f::BivariateFun)=space(f,1)⊗space(f,2)
domain(f::BivariateFun)=domain(f,1)*domain(f,2)

differentiate(u::BivariateFun,i::Integer,j::Integer)=j==0?u:differentiate(differentiate(u,i),i,j-1)
grad(u::BivariateFun)=[differentiate(u,1),differentiate(u,2)]
lap(u::BivariateFun)=differentiate(u,1,2)+differentiate(u,2,2)
Base.div{B<:BivariateFun}(u::Vector{B})=differentiate(u[1],1)+differentiate(u[2],2)
curl{B<:BivariateFun}(u::Vector{B})=differentiate(u[2],1)-differentiate(u[1],2)

Base.chop(f::MultivariateFun)=chop(f,10eps())
Base.eltype{T}(::MultivariateFun{T})=T
Base.eltype{T}(::Type{MultivariateFun{T}})=T
Base.eltype{MF<:MultivariateFun}(::Type{MF})=eltype(super(MF))


include("VectorFun.jl")
include("ProductDomain.jl")
include("TensorSpace.jl")
include("LowRankFun.jl")
include("ProductFun.jl")


arglength(f)=length(Base.uncompressed_ast(f.code.def).args[1])



## Convert between Fun and MultivariateFun
Fun(f::ProductFun)=Fun(fromtensor(coefficients(f)),space(f))
Fun(f::ProductFun,sp::TensorSpace)=Fun(ProductFun(f,sp))
Fun(f::LowRankFun)=Fun(ProductFun(f))


Fun(f::MultivariateFun,sp::FunctionSpace)=Fun(Fun(f),sp)

Fun(f,d1::Domain,d2::Domain)=Fun(f,d1*d2)
Fun{T<:Number,V<:Number}(f,d1::Vector{T},d2::Vector{V})=Fun(f,convert(Domain,d1),convert(Domain,d2))

coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])



points(f::BivariateFun,k...)=points(space(f),size(f,1),size(f,2),k...)


for OP in (:+,:-)
    @eval begin
        $OP(f::Fun,g::MultivariateFun)=$OP(ProductFun(f),g)
        $OP(f::MultivariateFun,g::Fun)=$OP(f,ProductFun(g))
    end
end


Base.sum{TS<:TensorSpace}(f::Fun{TS},k::Integer)=sum(ProductFun(f),k)
Base.sum{TS<:TensorSpace}(f::Fun{TS})=sum(ProductFun(f))


## kron

Base.kron(f::Fun,g::Fun)=Fun(LowRankFun([f],[g]))
Base.kron(f::Fun,g::Number)=kron(f,Fun(g))
Base.kron(f::Number,g::Fun)=kron(Fun(f),g)


