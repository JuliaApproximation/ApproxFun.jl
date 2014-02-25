export SingFun

type SingFun{T<:IFun} <: AbstractFun
    fun::T
    α::Float64
    β::Float64
end

Base.getindex(f::SingFun,x)=evaluate(f,x)

jacobiweight(α,β,x)=(1+x).^α.*(1-x).^β
jacobiweight(f::SingFun,x)=jacobiweight(f.α,f.β,tocanonical(f,x))





evaluate(f::SingFun,x)=f.fun[x].*jacobiweight(f,x)
values(f::SingFun)=values(f.fun).*jacobiweight(f,points(f))


for op in (:tocanonical,:fromcanonical)
    @eval begin
        ($op)(f::SingFun,x)=($op)(f.fun,x)
    end
end

for op in (:(Base.length),:points)
    @eval begin
        ($op)(f::SingFun)=($op)(f.fun)        
    end
end


pad(f::SingFun,n)=SingFun(pad(f.fun,n),f.α,f.β)