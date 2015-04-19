
export samplecdf,normalizedcumsum!


##bisection inverse


bisectioninv{S,T}(f::Fun{S,T},x::Real) = first(bisectioninv(f,[x]))


function bisectioninv{S,T}(f::Fun{S,T},x::Float64)
    d=domain(f)
    a = first(d);b = last(d)


    for k=1:47  #TODO: decide 47
        m=.5*(a+b)
        val = f[m]

            (val<= x) ? (a = m) : (b = m)
    end
    .5*(a+b)
end

bisectioninv{S,T}(f::Fun{S,T},x::Vector)=Float64[bisectioninv(f,xx) for xx in x]


## Clenshaw bisection

function chebbisectioninv(c::Vector{Float64},x::Float64)
    a = -1.;b = 1.


    for k=1:47  #TODO: decide 47
        m=.5*(a+b)
        val = clenshaw(c,m)

            (val<= x) ? (a = m) : (b = m)
    end
    .5*(a+b)
end


chebbisectioninv(c::Vector{Float64},xl::Vector{Float64})=(n=length(xl);chebbisectioninv(c,xl,ClenshawPlan(Float64,n)))
function chebbisectioninv(c::Vector{Float64},xl::Vector{Float64},plan::ClenshawPlan{Float64})
    n = length(xl)
    a = -ones(n)
    b = ones(n)


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
chebbisectioninv(c::Array{Float64,2},xl::Vector{Float64})=(n=length(xl);chebbisectioninv(c,xl,ClenshawPlan(Float64,n)))
function chebbisectioninv(c::Array{Float64,2},xl::Vector{Float64},plan::ClenshawPlan{Float64})
    @assert size(c)[2] == length(xl)

    n = length(xl)
    a = -ones(n)
    b = ones(n)


    for k=1:47  #TODO: decide 47
        m=.5*(a+b);
        vals = clenshaw(c,m,plan)

        for j = 1:n
            (vals[j] <= xl[j]) ? (a[j] = m[j]) : (b[j] = m[j])
        end
    end
    m=.5*(a+b)
end

for TYP in (:Vector,:Float64), SP in (:Chebyshev,:LineSpace)
    @eval bisectioninv(cf::Fun{$SP,Float64},x::$TYP)=fromcanonical(cf,chebbisectioninv(coefficients(cf),x))
end


##normalized cumsum

function normalizedcumsum(f::Fun)
    cf = cumsum(f)
    cf = cf/last(cf)

    cf
end

function subtract_zeroatleft!(f::Vector{Float64})
    for k=2:length(f)
        @inbounds f[1] += (-1.)^k.*f[k]
    end

    f
end

function subtract_zeroatleft!(f::Array{Float64,2})
    for k=2:size(f)[1],j=1:size(f)[2]
        @inbounds f[1,j] += (-1.)^k.*f[k,j]
    end

    f
end

function multiply_oneatright!(f::Vector{Float64})
    val=0.
    for k=1:length(f)
        val+=f[k]
    end

    val=1./val

    for k=1:length(f)
        @inbounds f[k] *= val
    end

    f
end

function multiply_oneatright!(f::Array{Float64,2})

    for j=1:size(f)[2]
        val=0.
        for k=1:size(f)[1]
            val+=f[k,j]
        end

        val=1./val

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
normalizedcumsum!(f::Vector{Float64})=chebnormalizedcumsum!(f)

## Sampling

sample(f::Fun,n::Integer)=samplecdf(normalizedcumsum(f),n)

samplecdf(cf::Fun,n::Integer)=bisectioninv(cf,rand(n))


sample(f::Fun)=sample(f,1)[1]
samplecdf(f::Fun)=samplecdf(f,1)[1]
samplecdf(v::Vector)=chebbisectioninv(v,rand())




##2D sample

sample{TS<:AbstractProductSpace}(f::Fun{TS},k::Integer)=sample(ProductFun(f),k)

function sample(f::LowRankFun{Chebyshev,Chebyshev,AbstractProductSpace{@compat(Tuple{Chebyshev,Chebyshev}),Float64,2},Float64,Float64},n::Integer)
    ry=sample(sum(f,1),n)
    fA=evaluate(f.A,ry)
    CB=coefficients(f.B)
    AB=CB*fA
    chebnormalizedcumsum!(AB)
    rx=chebbisectioninv(AB,rand(n))
  [fromcanonical(domain(f,1),rx) ry]
end



sample(f::MultivariateFun,n)=sample(LowRankFun(f),n)
sample(f::MultivariateFun)=sample(f,1)[1,:]




## Special spaces

# Rays may be schwartz at right endpoint so we project
function sample{SS<:JacobiWeight,DD<:Ray,TT}(f::Fun{MappedSpace{SS,DD,TT},Float64},n::Integer)
    if space(f).space.β == 0
        samplecdf(normalizedcumsum(f),n)
    else
        sample(Fun(x->f[x],MappedSpace(domain(f),JacobiWeight(space(f).space.α,1))),n)
    end
end

# Line/Ray have unbounded endpoints so we map
for TYP in (:Line,:Ray)
    @eval bisectioninv{SS,DD<:$TYP,TT}(f::Fun{MappedSpace{SS,DD,TT},Float64},x::Vector)=fromcanonical(f,bisectioninv(Fun(f.coefficients,space(f).space),x))
end



function sample{SS}(f::LowRankFun{LineSpace{SS},LineSpace{SS},AbstractProductSpace{@compat(Tuple{LineSpace{SS},LineSpace{SS}}),Float64,2},Float64,Float64},n::Integer)
    cf=normalizedcumsum(sum(f,1))
    CB=coefficients(map(cumsum,f.B))

    ry=samplecdf(cf,n)
    fA=evaluate(f.A,ry)
    CBfA=CB*fA  #cumsums at points
    multiply_oneatright!(CBfA)

    rx=fromcanonical(first(f.B),chebbisectioninv(CBfA,rand(n)))

    [rx ry]
end


