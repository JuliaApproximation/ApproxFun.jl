abstract type MultivariateFun{T,N} end
const BivariateFun{T} = MultivariateFun{T,2}

export grad, lap, curl

#implements coefficients/values/evaluate
space{T,N}(f::MultivariateFun{T,N})=mapreduce(k->space(f,k),⊗,1:N)
domain{T,N}(f::MultivariateFun{T,N})=mapreduce(k->domain(f,k),*,1:N)

domain(f::MultivariateFun,k::Integer)=domain(space(f,k))

differentiate(u::BivariateFun,i::Integer,j::Integer) =
    j==0?u:differentiate(differentiate(u,i),i,j-1)
grad(u::BivariateFun) = [differentiate(u,1),differentiate(u,2)]
lap(u::BivariateFun) = differentiate(u,1,2)+differentiate(u,2,2)
Base.div{B<:BivariateFun}(u::Vector{B}) =
    differentiate(u[1],1)+differentiate(u[2],2)
curl{B<:BivariateFun}(u::Vector{B}) = differentiate(u[2],1)-differentiate(u[1],2)

Base.chop(f::MultivariateFun) = chop(f,10eps())
Base.eltype{T}(::MultivariateFun{T}) = T
Base.eltype{T,N}(::Type{MultivariateFun{T,N}}) = T
Base.eltype{MF<:MultivariateFun}(::Type{MF}) = eltype(supertype(MF))

include("VectorFun.jl")
include("TensorSpace.jl")
include("LowRankFun.jl")
include("ProductFun.jl")


arglength(f)=length(Base.uncompressed_ast(f.code.def).args[1])



## Convert between Fun and MultivariateFun
# need to chop extra zeros
Fun(f::ProductFun) =
    Fun(space(f),chop!(fromtensor(space(f),coefficients(f)),0))
Fun(f::ProductFun,sp::TensorSpace) = Fun(ProductFun(f,sp))
Fun(f::LowRankFun) = Fun(ProductFun(f))


Fun(f::MultivariateFun,sp::Space) = Fun(Fun(f),sp)

Fun(f,d1::Domain,d2::Domain) = Fun(f,d1*d2)
Fun{T<:Number,V<:Number}(f,d1::Vector{T},d2::Vector{V}) = Fun(f,convert(Domain,d1),convert(Domain,d2))

coefficients(f::BivariateFun,sp::TensorSpace)=coefficients(f,sp[1],sp[2])



points(f::BivariateFun,k...)=points(space(f),size(f,1),size(f,2),k...)


function *(vx::LowRankFun,u0::ProductFun)
    ret=zeros(space(u0))
    for k=1:length(vx.A)
        a,b=vx.A[k],vx.B[k]
        ret+=((b*((a*u0).')).')
    end
    ret
end

*(a::ProductFun,b::LowRankFun)=b*a
*(a::MultivariateFun,b::MultivariateFun)=LowRankFun(a)*ProductFun(b)

for OP in (:+,:-,:*,:/)
    @eval begin
        $OP(f::Fun,g::MultivariateFun)=$OP(ProductFun(f),g)
        $OP(f::MultivariateFun,g::Fun)=$OP(f,ProductFun(g))
    end
end


Base.sum{TS<:TensorSpace}(f::Fun{TS},k::Integer)=sum(ProductFun(f),k)
Base.sum{TS<:TensorSpace}(f::Fun{TS})=sum(ProductFun(f))


## kron
# TODO: generalize
function Base.kron(f::Fun,g::Fun)
    sp=space(f)⊗space(g)
    it=tensorizer(sp)
    N=ncoefficients(f);M=ncoefficients(g)
    cfs=Array{promote_type(eltype(f),eltype(g))}(0)
    for (k,j) in it
        # Tensor product is N x M, so if we are outside
        # the (N+M)th diagonal we have no more entries
        if k+j > N+M
            break
        elseif k ≤ N && j ≤ M
            push!(cfs,f.coefficients[k]*g.coefficients[j])
        else
            push!(cfs,0)
        end
    end
    Fun(sp,cfs)
end
Base.kron(f::Fun,g::Number) = kron(f,Fun(g))
Base.kron(f::Number,g::Fun) = kron(Fun(f),g)
