


ToeplitzOperator{DD}(f::Fun{Laurent{DD}})=ToeplitzOperator(f.coefficients[2:2:end],
                                                    f.coefficients[1:2:end])



##Taylor

Multiplication{DD}(f::Fun{Taylor{DD}},sp::Taylor{DD}) =
    MultiplicationWrapper(f,SpaceOperator(ToeplitzOperator(f.coefficients[2:end],
                                                           [f.coefficients[1]]),
                                          sp,sp))


## Evaluation

getindex{DD}(T::ConcreteEvaluation{Taylor{DD},Complex{Float64},Int,Complex{Float64}},cols::Range) =
    mappoint(domain(T),Circle(),T.x).^(cols-1)


## Multiplication

Multiplication{DD}(f::Fun{Laurent{DD}},sp::Laurent{DD}) = ConcreteMultiplication(eltype(f),f,sp)

function laurent_getindex{T}(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer)
    # switch to double-infinite indices
    k=iseven(k)?-k÷2:(k-1)÷2
    j=iseven(j)?-j÷2:(j-1)÷2

    if 0<k-j≤length(negative)
        negative[k-j]
    elseif 0≤j-k≤length(nonnegative)-1
        nonnegative[j-k+1]
    else
        zero(T)
    end
end

rangespace{DD}(T::ConcreteMultiplication{Laurent{DD},Laurent{DD}}) = domainspace(T)
function getindex{DD}(T::ConcreteMultiplication{Laurent{DD},Laurent{DD}},k::Integer,j::Integer)
    isempty(T.f.coefficients) && return zero(eltype(T))
    laurent_getindex(T.f.coefficients[3:2:end],T.f.coefficients[[1;2:2:end]],k,j)
end

function bandinds{DD}(T::ConcreteMultiplication{Laurent{DD},Laurent{DD}})
    bbi = blockbandinds(T)
    (2bbi[1]-1,2bbi[2]+1)
end

function blockbandinds{DD}(T::ConcreteMultiplication{Laurent{DD},Laurent{DD}})
    m = ncoefficients(T.f)÷2
    (-m,m)
end


Multiplication{DD}(f::Fun{Fourier{DD}},sp::Laurent{DD}) = Multiplication(Fun(f,sp),sp)
Multiplication{DD}(f::Fun{Laurent{DD}},sp::Fourier{DD}) = Multiplication(Fun(f,sp),sp)

# override SumSpace default
coefficienttimes{DD}(f::Fun{Laurent{DD}},g::Fun{Laurent{DD}}) = Multiplication(f,space(g))*g




## Derivative

# override map definition
Derivative{s,DD<:Circle}(S::Hardy{s,DD},k::Integer) = ConcreteDerivative(S,k)
Derivative{s,DD<:PeriodicInterval}(S::Hardy{s,DD},k::Integer) = ConcreteDerivative(S,k)
Derivative{DD<:Circle}(S::Laurent{DD},k::Integer) =
    DerivativeWrapper(InterlaceOperator(Diagonal([map(s->Derivative(s,k),S.spaces)...]),SumSpace),k)

bandinds{s,DD<:PeriodicInterval}(D::ConcreteDerivative{Hardy{s,DD}})=(0,0)
bandinds{s,DD<:Circle}(D::ConcreteDerivative{Hardy{s,DD}})=s?(0,D.order):(-D.order,0)

rangespace{S<:Hardy}(D::ConcreteDerivative{S})=D.space

function taylor_derivative_getindex(d::PeriodicInterval,m,k::Integer,j::Integer)
    C=2/(d.b-d.a)*π*im
    k==j?(C*(k-1))^m:zero(C)
end

function hardyfalse_derivative_getindex(d::PeriodicInterval,m,k::Integer,j::Integer)
    C=2/(d.b-d.a)*π*im
    k==j?(-C*k)^m:zero(C)
end



function taylor_derivative_getindex(d::Circle,m,k::Integer,j::Integer)
    C=d.radius^(-m)

    if j==k+m
        D=k
        for j=k+1:k+m-1
          D*=j
        end
        C*D
    else
        zero(C)
    end
end

function hardyfalse_derivative_getindex(d::Circle,m,k::Integer,j::Integer)
    C=(-d.radius)^(-m)

    if j == k-m
        D=k-m
        for j=k-m+1:k-1
          D*=j
        end
        C*D
    else
        zero(C)
    end
end


getindex{DD,OT,T}(D::ConcreteDerivative{Taylor{DD},OT,T},k::Integer,j::Integer)  =
    T(taylor_derivative_getindex(domain(D),D.order,k,j))
getindex{DD,OT,T}(D::ConcreteDerivative{Hardy{false,DD},OT,T},k::Integer,j::Integer) =
    T(hardyfalse_derivative_getindex(domain(D),D.order,k,j))




## Integral



# Integral{D<:Circle}(S::Taylor{D},m)=Integral{Taylor,typeof(m),Complex{Float64}}(S,m)
#
# function Integral{D<:PeriodicInterval}(S::Hardy{false,D},m)
#     Integral{Hardy{false,D},typeof(m),Complex{Float64}}(S,m)
# end


Integral{s,DD<:Circle}(S::Hardy{s,DD},k::Integer)=ConcreteIntegral(S,k)

bandinds{DD<:Circle}(D::ConcreteIntegral{Taylor{DD}})=(-D.order,0)
rangespace{s,DD<:Circle}(Q::ConcreteIntegral{Hardy{s,DD}})=Q.space

function getindex{DD<:Circle}(D::ConcreteIntegral{Taylor{DD}},k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    T=eltype(D)
    C=T(d.radius^m)

    if j==k-m
        D=k-m
        for j=k-m+1:k-1
          D*=j
        end
        C/D
    else
        zero(T)
    end
end


function Integral{T,DD<:Circle}(S::SubSpace{Hardy{false,DD},UnitCount{Int64},T,DD,1},k::Integer)
    if first(S.indexes) == k+1
        ConcreteIntegral(S,k)
    else
        @assert first(S.index) > k+1
        S2=S.space|(k+1:∞)
        IntegralWrapper(ConcreteIntegral(S2,k)*Converion(S,S2),k)
    end
end

bandinds{T,DD<:Circle}(D::ConcreteIntegral{SubSpace{Hardy{false,DD},UnitCount{Int64},T,DD,1}}) =
    (0,0)

rangespace{T,DD<:Circle}(D::ConcreteIntegral{SubSpace{Hardy{false,DD},UnitCount{Int64},T,DD,1}}) =
    D.space.space

function getindex{T,OT,TT,DD<:Circle}(D::ConcreteIntegral{SubSpace{Hardy{false,DD},UnitCount{Int64},T,DD,1},OT,TT},
                                      k::Integer,j::Integer)
    d=domain(D)
    m=D.order

    C=TT((-d.radius)^m)

    if k==j
        D=k
        for j=k+1:k+m-1
          D*=j
        end
        C/D
    else
        zero(TT)
    end
end



bandinds{DD<:PeriodicInterval}(D::ConcreteIntegral{Hardy{false,DD}}) = (0,0)
rangespace{DD<:PeriodicInterval}(D::ConcreteIntegral{Taylor{DD}}) = D.space


function getindex{DD<:PeriodicInterval}(D::ConcreteIntegral{Hardy{false,DD}},k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    T=eltype(D)
    C=2/(d.b-d.a)*π*im
    if k==j
        T((-C*k)^(-m))
    else
        zero(T)
    end
end



bandinds{T,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{Taylor{DD},UnitCount{Int64},T,DD,1}}) =
    (0,0)
rangespace{T,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{Taylor{DD},UnitCount{Int64},T,DD,1}}) =
    D.space

function getindex{T,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{Taylor{DD},UnitCount{Int64},T,DD,1}},
                                          k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    TT=eltype(D)
    C=2/(d.b-d.a)*π*im
    if k==j
        TT((C*(k+n-1))^(-m))
    else
        zero(TT)
    end
end



## Definite integral

for SP in (:Taylor,:(Hardy{false}),:Laurent)
    @eval begin
        DefiniteIntegral{D}(S::$SP{D}) =
            ConcreteDefiniteIntegral{typeof(S),promote_type(eltype(S),eltype(D))}(S)
        DefiniteLineIntegral{D}(S::$SP{D}) =
            ConcreteDefiniteLineIntegral{typeof(S),real(promote_type(eltype(S),eltype(D)))}(S)
    end
end

getindex{T,D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Taylor{D},T},k::Integer) =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex{T,D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Hardy{false,D},T},k::Integer) =
    zero(T)

getindex{T,D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Laurent{D},T},k::Integer) =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex{T,D<:Circle}(Σ::ConcreteDefiniteIntegral{Taylor{D},T},k::Integer) =
    zero(T)

getindex{T,D<:Circle}(Σ::ConcreteDefiniteIntegral{Hardy{false,D},T},k::Integer) =
    k == 1? T(complexlength(domain(Σ))) :zero(T)

getindex{T,D<:Circle}(Σ::ConcreteDefiniteIntegral{Laurent{D},T},k::Integer) =
    k == 2? T(complexlength(domain(Σ))) :zero(T)

getindex{T,D}(Σ::ConcreteDefiniteLineIntegral{Taylor{D},T},k::Integer) =
    k == 1? T(arclength(domain(Σ))) : zero(T)

getindex{T,D}(Σ::ConcreteDefiniteLineIntegral{Hardy{false,D},T},k::Integer) =
    zero(T)

getindex{T,D}(Σ::ConcreteDefiniteLineIntegral{Laurent{D},T},k::Integer) =
    k == 1? T(arclength(domain(Σ))) : zero(T)

bandinds{D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Taylor{D}}) = 0,0
bandinds{D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Hardy{false,D}}) = 0,0
bandinds{D<:PeriodicInterval}(Σ::ConcreteDefiniteIntegral{Laurent{D}}) = 0,0
bandinds{D<:Circle}(Σ::ConcreteDefiniteIntegral{Taylor{D}}) = 0,0
bandinds{D<:Circle}(Σ::ConcreteDefiniteIntegral{Hardy{false,D}}) = 0,0
bandinds{D<:Circle}(Σ::ConcreteDefiniteIntegral{Laurent{D}}) = 0,1
bandinds{D}(Σ::ConcreteDefiniteLineIntegral{Taylor{D}}) = 0,0
bandinds{D}(Σ::ConcreteDefiniteLineIntegral{Hardy{false,D}}) = 0,0
bandinds{D}(Σ::ConcreteDefiniteLineIntegral{Laurent{D}}) = 0,0


## reverse orientation

conversion_type{DD<:Circle}(A::Laurent{DD},B::Laurent{DD})=domain(A).orientation?A:B
function Conversion{DD}(A::Laurent{DD},B::Laurent{DD})
    @assert domain(A) == reverse(domain(B))
    ConversionWrapper(SpaceOperator(
        InterlaceOperator(Diagonal([eye(1),PermutationOperator([2,1])]))
    ,A,B))
end
