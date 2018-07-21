
export samplecdf,normalizedcumsum!


##bisection inverse


bisectioninv(f::Fun{S,T},x::Real;opts...) where {S,T} = first(bisectioninv(f,[x];opts...))


function bisectioninv(f::Fun{S,T},x::Float64;numits::Int=47) where {S,T}
    d=domain(f)
    a = first(d);b = last(d)


    for k=1:numits  #TODO: decide 47
        m=.5*(a+b)
        val = f(m)

            (val<= x) ? (a = m) : (b = m)
    end
    .5*(a+b)
end

bisectioninv(f::Fun{S,T},x::AbstractVector;opts...) where {S,T} = Float64[bisectioninv(f,xx;opts...) for xx in x]


## Clenshaw bisection

function chebbisectioninv(c::AbstractVector{Float64},x::Float64;numits::Int=47)
    a = -1.;b = 1.

    C=Chebyshev()
    for k=1:numits  #TODO: decide 47
        m=.5*(a+b)
        val = clenshaw(C,c,m)

            (val<= x) ? (a = m) : (b = m)
    end
    .5*(a+b)
end


chebbisectioninv(c::AbstractVector{Float64},xl::AbstractVector{Float64}) =
    (n=length(xl);chebbisectioninv(c,xl,ClenshawPlan(Float64,Chebyshev(),length(c),n)))
function chebbisectioninv(c::AbstractVector{Float64},xl::AbstractVector{Float64},plan::ClenshawPlan{Chebyshev{D,R},Float64}) where {D<:Domain,R}
    n = length(xl)
    a = -fill(1.0,n)
    b = fill(1.0,n)


    for k=1:47  #TODO: decide 47
        m=.5*(a+b);
        vals = clenshaw(c,m,plan)

        for j = 1:n
            (vals[j] <= xl[j]) ? (a[j] = m[j]) : (b[j] = m[j])
        end
    end
    m=.5*(a+b)
end


#here, xl is vector w/ length == #cols of c
chebbisectioninv(c::AbstractMatrix{Float64},xl::AbstractVector{Float64}) =
    (n=length(xl);chebbisectioninv(c,xl,ClenshawPlan(Float64,Chebyshev(),size(c,1),n)))
function chebbisectioninv(c::AbstractMatrix{Float64},xl::AbstractVector{Float64},plan::ClenshawPlan{Chebyshev{D,R},Float64}) where {D<:Domain,R}
    @assert size(c)[2] == length(xl)

    n = length(xl)
    a = -fill(1.0,n)
    b = fill(1.0,n)


    for k=1:47  #TODO: decide 47
        m=.5*(a+b);
        vals = clenshaw(c,m,plan)

        for j = 1:n
            (vals[j] <= xl[j]) ? (a[j] = m[j]) : (b[j] = m[j])
        end
    end
    m=.5*(a+b)
end

for TYP in (:Vector,:Float64)
    @eval begin
        bisectioninv(cf::Fun{SP,Float64},x::$TYP;opts...) where {SP<:Chebyshev} =
            fromcanonical.(Ref(space(cf)),chebbisectioninv(coefficients(cf),x;opts...))
#        bisectioninv{SP<:LineSpace}(cf::Fun{SP,Float64},x::$TYP;opts...)=fromcanonical(cf,chebbisectioninv(coefficients(cf),x;opts...))
    end
end


##normalized cumsum

function normalizedcumsum(f::Fun)
    cf = cumsum(f)
    cf = cf/last(cf)

    cf
end

function subtract_zeroatleft!(f::AbstractVector{Float64})
    for k=2:length(f)
        @inbounds f[1] += (-1.)^k.*f[k]
    end

    f
end

function subtract_zeroatleft!(f::AbstractMatrix{Float64})
    for k=2:size(f)[1],j=1:size(f)[2]
        @inbounds f[1,j] += (-1.)^k.*f[k,j]
    end

    f
end

function multiply_oneatright!(f::AbstractVector{Float64})
    val=0.0
    for k=1:length(f)
        val += f[k]
    end

    val = 1/val

    for k=1:length(f)
        @inbounds f[k] *= val
    end

    f
end

function multiply_oneatright!(f::AbstractMatrix{Float64})

    for j=1:size(f)[2]
        val=0.0
        for k=1:size(f)[1]
            val += f[k,j]
        end

        val = 1/val

        for k=1:size(f)[1]
            @inbounds f[k,j] *= val
        end
    end

    f
end

function chebnormalizedcumsum!(f)
    ultraconversion!(f)
    ultraint!(f)
    subtract_zeroatleft!(f)
    multiply_oneatright!(f)
end

# For RandomMatrices compatibility
normalizedcumsum!(f::AbstractVector{Float64})=chebnormalizedcumsum!(f)

## Sampling

sample(f::Fun,n::Integer) = samplecdf(normalizedcumsum(f),n)

samplecdf(cf::Fun,n::Integer) = bisectioninv(cf,rand(n))


sample(f::Fun) = sample(f,1)[1]
samplecdf(f::Fun) = samplecdf(f,1)[1]
samplecdf(v::AbstractVector) = chebbisectioninv(v,rand())




##2D sample

sample(f::Fun{TS},k::Integer) where {TS<:AbstractProductSpace}  =sample(ProductFun(f),k)

function sample(f::LowRankFun,n::Integer)
    rx=sample(sum(f,2),n)
    fA=evaluate.(f,rx,:)
    ry=map(sample,fA)
    [rx ry]
end

function sample(f::LowRankFun{C,C,TensorSpace{Tuple{C,C},DD,RR},Float64},n::Integer) where {C<:Chebyshev,DD<:BivariateDomain,RR<:Real}
    ry=sample(sum(f,1),n)
    fA=evaluate.(f.A,ry')
    CB=coefficientmatrix(f.B)
    AB=CB*fA
    chebnormalizedcumsum!(AB)
    rx=chebbisectioninv(AB,rand(n))
    [fromcanonical.(Ref(domain(f,1)),rx) ry]
end



sample(f::ProductFun,n)=sample(LowRankFun(f),n)
sample(f::MultivariateFun)=sample(f,1)[1,:]




## Special spaces

# Rays may be schwartz at right endpoint so we project
function sample(f::Fun{JacobiWeight{SS,DD,RR,TT},Float64}, n::Integer) where {SS,DD<:Ray,RR,TT}
    if space(f).α == 0
        samplecdf(normalizedcumsum(f),n)
    else
        sample(default_Fun(f,JacobiWeight(1,space(f).α,domain(f))),n)
    end
end
