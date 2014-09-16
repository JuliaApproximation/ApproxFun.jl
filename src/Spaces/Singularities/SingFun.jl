export SingFun



type SingFun{T<:Fun} <: AbstractFun
    fun::T
    α::Float64
    β::Float64
end

SingFun(f,a::Integer,b::Integer)=SingFun(f,1.a,1.b)

SingFun(f::Fun)=SingFun(f,0.,0.)

Base.getindex(f::SingFun,x)=evaluate(f,x)

jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
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
    end
end

for op in (:*,:.*)
    @eval begin
        ($op)(c::Number,f::SingFun)=SingFun(($op)(c,f.fun),f.α,f.β)        
    end
end

for op in (:/,:./)
    @eval begin
        ($op)(c::Number,f::SingFun)=SingFun(($op)(c,f.fun),-f.α,-f.β)        
    end
end

.*(f::SingFun,g::SingFun)=SingFun(f.fun.*g.fun,f.α+g.α,f.β+g.β)
./(f::SingFun,g::SingFun)=SingFun(f.fun./g.fun,f.α-g.α,f.β-g.β)

for op in (:./,:.*)
    @eval begin
        ($op)(f::SingFun,g::Fun)=($op)(f,SingFun(g))
        ($op)(f::Fun,g::SingFun)=($op)(SingFun(f),g)
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

## transform alpha, beta
# We assume that the user knows whether this is possible
#TODO: We also assume it is applied to "Schwartz data", but this can be fixed


function increase_jacobi_parameter(s,f::SingFun)
    if s == -1
        SingFun(divide_singularity(s,f.fun),f.α+1,f.β)
    elseif s == 1
        SingFun(divide_singularity(s,f.fun),f.α,f.β+1)    
    end
end

increase_jacobi_parameter(f::SingFun)=SingFun(divide_singularity(f.fun),f.α+1,f.β+1)    




## Calculus

function Base.sum(f::SingFun)
    ##TODO: generalize

    if f.α==f.β==.5
        fromcanonicalD(f,0.)*coefficients(f.fun,UltrasphericalSpace{1}(domain(f)))[1]*π/2
    elseif f.α==f.β==0.
        sum(f.fun)
    elseif f.α==f.β==-.5
        fromcanonicalD(f,0.)*π*f.fun.coefficients[1]
    elseif f.α<0. && f.β<0.
        #TODO: should be < -1.
        sum(increase_jacobi_parameter(f))
    elseif f.α < 0
        sum(increase_jacobi_parameter(-1,f))
    elseif  f.β < 0
        sum(increase_jacobi_parameter(+1,f))    
    else
        error("sum not implemented for all Jacobi parameters")
    end
end



function integrate(f::SingFun)    
    if f.α==f.β==0.
        SingFun(integrate(f.fun),0.,0.)
    else
        d=domain(f)
            
        @assert d==Interval()
        @assert f.α==f.β==-.5
    
        SingFun(Fun(ultraiconversion(-f.fun.coefficients[2:end]./[1:length(f)-1])),.5,.5)
    end
end
