
export Chebyshev


immutable Chebyshev{D<:Domain} <: PolynomialSpace{D}
    domain::D
    Chebyshev(d) = new(d)
    Chebyshev() = new(Interval())
end

Chebyshev() = Chebyshev{Interval{Float64}}()
Chebyshev(d::Domain) = Chebyshev{typeof(d)}(d)
Chebyshev(d::Vector) = Chebyshev(Domain(d))


Space(d::Interval) = Chebyshev(d)


setdomain(S::Chebyshev,d::Domain) = Chebyshev(d)

Base.ones{T<:Number}(::Type{T},S::Chebyshev) = Fun(ones(T,1),S)
Base.ones(S::Chebyshev) = Fun(ones(1),S)

Base.first{D}(f::Fun{Chebyshev{D}}) = foldr(-,coefficients(f))
Base.last{D}(f::Fun{Chebyshev{D}}) = reduce(+,coefficients(f))
identity_fun(d::Chebyshev) = identity_fun(domain(d))

spacescompatible(a::Chebyshev,b::Chebyshev) = domainscompatible(a,b)
hasfasttransform(::Chebyshev) = true


function coefficients(g::Vector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function coefficients(g::Vector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::Vector,plan) = chebyshevtransform(vals,plan)
itransform(::Chebyshev,cfs::Vector,plan) = ichebyshevtransform(cfs,plan)
plan_transform(::Chebyshev,vals::Vector) = plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::Vector) = plan_ichebyshevtransform(cfs)

## Evaluation

clenshaw(sp::Chebyshev,c::AbstractVector,x::AbstractArray) =
    clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))

function clenshaw(::Chebyshev,c::AbstractVector,x)
    N,T = length(c),promote_type(eltype(c),typeof(x))
    if isempty(c)
        return zero(x)
    end

    x = 2x
    bk1,bk2 = zero(T),zero(T)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd(x,bk1,c[k]-bk2)
    end

    muladd(x/2,bk1,c[1]-bk2)
end

clenshaw{S<:Chebyshev,V}(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V})=
    clenshaw(c,collect(x),plan)

#TODO: This modifies x, which is not threadsafe
function clenshaw{S<:Chebyshev,V}(c::AbstractVector,x::Vector,plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(V,n)
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    @inbounds for i = 1:n
        x[i] = 2x[i]
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck = c[k]
        for i = 1:n
            bk[i] = muladd(x[i],bk1[i],ck-bk2[i])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ck = c[1]
    @inbounds for i = 1:n
        x[i] = x[i]/2
        bk[i] = muladd(x[i],bk1[i],ck-bk2[i])
    end

    bk
end

function clenshaw{S<:Chebyshev,T<:Number}(c::AbstractMatrix{T},x::T,plan::ClenshawPlan{S,T})
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    m,n=size(c) # m is # of coefficients, n is # of funs

    @inbounds for i = 1:n
        bk1[i] = zero(T)
        bk2[i] = zero(T)
    end
    x = 2x

    @inbounds for k=m:-1:2
        for j=1:n
            ck = c[k,j]
            bk[j] = muladd(x,bk1[j],ck - bk2[j])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    x = x/2
    @inbounds for i = 1:n
        ce = c[1,i]
        bk[i] = muladd(x,bk1[i],ce - bk2[i])
    end

    bk
end

function clenshaw{S<:Chebyshev,T<:Number}(c::AbstractMatrix{T},x::AbstractVector{T},plan::ClenshawPlan{S,T})
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    m,n=size(c) # m is # of coefficients, n is # of funs

    @inbounds for i = 1:n
        x[i] = 2x[i]
        bk1[i] = zero(T)
        bk2[i] = zero(T)
    end

    @inbounds for k=m:-1:2
        for j=1:n
            ck = c[k,j]
            bk[j] = muladd(x[j],bk1[j],ck - bk2[j])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end


    @inbounds for i = 1:n
        x[i] = x[i]/2
        ce = c[1,i]
        bk[i] = muladd(x[i],bk1[i],ce - bk2[i])
    end

    bk
end

# overwrite x

function clenshaw!{S<:Chebyshev,V}(c::Vector,x::Vector,plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)

    if isempty(c)
        fill!(x,0)
        return x
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    @inbounds for i = 1:n
        x[i] = 2x[i]
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck = c[k]
        for i = 1:n
            bk[i] = muladd(x[i],bk1[i],ck - bk2[i])
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ce = c[1]
    @inbounds for i = 1:n
        x[i] = x[i]/2
        x[i] = muladd(x[i],bk1[i],ce-bk2[i])
    end

    x
end

## Calculus


# diff T -> U, then convert U -> T
integrate{D<:Interval}(f::Fun{Chebyshev{D}}) =
    Fun(fromcanonicalD(f,0)*ultraint!(ultraconversion(f.coefficients)),f.space)
differentiate{D<:Interval}(f::Fun{Chebyshev{D}}) =
    Fun(1/fromcanonicalD(f,0)*ultraiconversion(ultradiff(f.coefficients)),f.space)

## identity_fun





## Multivariate


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op){S<:Chebyshev,V<:Chebyshev}(f::ProductFun{S,V}) =
        ProductFun(chebyshevtransform($op(values(f))),space(f))
end



reverseorientation{C<:Chebyshev}(f::Fun{C}) =
    Fun(alternatesign!(copy(f.coefficients)),Chebyshev(reverse(domain(f))))


include("ChebyshevOperators.jl")
