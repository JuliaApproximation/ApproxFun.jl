export SingFun

type SingFun{T<:IFun} <: AbstractFun
    fun::T
    α::Float64
    β::Float64
end

SingFun(f::IFun)=SingFun(f,0.,0.)

Base.getindex(f::SingFun,x)=evaluate(f,x)

jacobiweight(α,β,x)=(1+x).^α.*(1-x).^β
jacobiweight(f::SingFun,x)=jacobiweight(f.α,f.β,tocanonical(f,x))





evaluate(f::SingFun,x)=f.fun[x].*jacobiweight(f,x)
values(f::SingFun)=values(f.fun).*jacobiweight(f,points(f))


for op in (:tocanonical,:fromcanonical,:fromcanonicalD,:tocanonicalD)
    @eval begin
        ($op)(f::SingFun,x)=($op)(f.fun,x)
    end
end

for op in (:(Base.length),:points,:domain)
    @eval begin
        ($op)(f::SingFun)=($op)(f.fun)        
    end
end

for op in (:*,:.*,:/,:./)
    @eval begin
        ($op)(f::SingFun,c::Number)=SingFun(($op)(f.fun,c),f.α,f.β)
        ($op)(c::Number,f::SingFun)=SingFun(($op)(c,f.fun),f.α,f.β)        
    end
end

.*(f::SingFun,g::SingFun)=SingFun(f.fun.*g.fun,f.α+g.α,f.β+g.β)
./(f::SingFun,g::SingFun)=SingFun(f.fun./g.fun,f.α-g.α,f.β-g.β)

for op in (:./,:.*)
    @eval begin
        ($op)(f::SingFun,g::IFun)=($op)(f,SingFun(g))
        ($op)(f::IFun,g::SingFun)=($op)(SingFun(f),g)
    end
end

for op in (:+,:-)
    @eval begin
        function ($op)(f::SingFun,g::SingFun)
            @assert f.α==g.α
            @assert f.β==g.β
            SingFun(($op)(f.fun,g.fun),f.α,f.β)
        end
    end
end



pad(f::SingFun,n)=SingFun(pad(f.fun,n),f.α,f.β)


## Calculus

function Base.sum(f::SingFun)
    ##TODO: generalize

    @assert f.α==.5
    @assert f.β==.5
    fromcanonicalD(f,0.)*coefficients(f.fun,1)[1]*π/2
end

