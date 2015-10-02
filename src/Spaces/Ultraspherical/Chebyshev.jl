
export Chebyshev


typealias Chebyshev{D<:Domain} Ultraspherical{0,D}


Space(d::Interval)=Chebyshev(d)
canonicalspace(S::Ultraspherical)=Chebyshev(domain(S))

function coefficients(g::Vector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function coefficients(g::Vector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::Vector,plan)=chebyshevtransform(vals,plan)
itransform(::Chebyshev,cfs::Vector,plan)=ichebyshevtransform(cfs,plan)
plan_transform(::Chebyshev,vals::Vector)=plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::Vector)=plan_ichebyshevtransform(cfs)

## Evaluation

function clenshaw{U<:Number,V<:Number}(::Chebyshev,c::AbstractVector{U},x::V)
    N,T = length(c),promote_type(U,V)
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

function clenshaw{S<:Chebyshev,T<:Number,U<:Number,V<:Number}(c::AbstractVector{T},x::AbstractVector{U},plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(x)
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

function clenshaw!{S<:Chebyshev,T<:Number,U<:Number,V<:Number}(c::Vector{T},x::Vector{U},plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)

    if isempty(c)
        for k=1:n
            x[k]=zero(U)
        end
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
integrate{C<:Chebyshev}(f::Fun{C})=Fun(chebyshevintegrate(domain(f),f.coefficients),f.space)
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint!(ultraconversion(cfs))


differentiate{C<:Chebyshev}(f::Fun{C})=Fun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(Fun(x->tocanonicalD(d,x),d).*Fun(differentiate(Fun(cfs)),d)).coefficients


## identity_fun

identity_fun(d::Chebyshev)=identity_fun(domain(d))


## Piecewise union

# union_rule dictates how to create a space that both spaces can be converted to
# in this case, it means
function union_rule{S1<:Tuple{Vararg{Ultraspherical}},
                    S2<:Tuple{Vararg{Ultraspherical}}}(s1::PiecewiseSpace{S1},s2::PiecewiseSpace{S2})
    PiecewiseSpace(map(Chebyshev,merge(domain(s1),domain(s2)).domains))
end

function union_rule{S1<:Tuple{Vararg{Ultraspherical}}}(s1::PiecewiseSpace{S1},s2::Ultraspherical)
    PiecewiseSpace(map(Chebyshev,merge(domain(s1),domain(s2)).domains))
end



## Multivariate


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op){S<:Chebyshev,V<:Chebyshev}(f::ProductFun{S,V})=ProductFun(chebyshevtransform($op(values(f))),space(f))
end



reverseorientation{C<:Chebyshev}(f::Fun{C})=Fun(alternatesign!(copy(f.coefficients)),Chebyshev(reverse(domain(f))))
