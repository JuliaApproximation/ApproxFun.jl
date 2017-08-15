


ToeplitzOperator(f::Fun{Laurent{DD,RR}}) where {DD,RR}=ToeplitzOperator(f.coefficients[2:2:end],
                                                    f.coefficients[1:2:end])



##Taylor

Multiplication(f::Fun{Taylor{DD,RR}},sp::Taylor{DD,RR}) where {DD,RR} =
    MultiplicationWrapper(f,SpaceOperator(ToeplitzOperator(f.coefficients[2:end],
                                                           [f.coefficients[1]]),
                                          sp,sp))


## Evaluation

getindex(T::ConcreteEvaluation{Taylor{DD,RR},Complex{Float64},Int,Complex{Float64}},cols::Range) where {DD,RR} =
    mappoint(domain(T),Circle(),T.x).^(cols-1)


## Multiplication

Multiplication(f::Fun{Laurent{DD,RR}},sp::Laurent{DD,RR}) where {DD,RR} = ConcreteMultiplication(eltype(f),f,sp)

function laurent_getindex(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer) where T
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

rangespace(T::ConcreteMultiplication{Laurent{DD,RR},Laurent{DD,RR}}) where {DD,RR} = domainspace(T)
function getindex(T::ConcreteMultiplication{Laurent{DD,RR},Laurent{DD,RR}},k::Integer,j::Integer) where {DD,RR}
    isempty(T.f.coefficients) && return zero(eltype(T))
    laurent_getindex(T.f.coefficients[3:2:end],T.f.coefficients[[1;2:2:end]],k,j)
end

function bandinds(T::ConcreteMultiplication{Laurent{DD,RR},Laurent{DD,RR}}) where {DD,RR}
    bbi = blockbandinds(T)
    (2bbi[1]-1,2bbi[2]+1)
end

function blockbandinds(T::ConcreteMultiplication{Laurent{DD,RR},Laurent{DD,RR}}) where {DD,RR}
    m = ncoefficients(T.f)÷2
    (-m,m)
end


Multiplication(f::Fun{Fourier{DD,R1}},sp::Laurent{DD,R2}) where {DD,R1,R2} = Multiplication(Fun(f,sp),sp)
Multiplication(f::Fun{Laurent{DD,R1}},sp::Fourier{DD,R2}) where {DD,R1,R2} = Multiplication(Fun(f,sp),sp)

# override SumSpace default
coefficienttimes(f::Fun{Laurent{DD,RR}},g::Fun{Laurent{DD,RR}}) where {DD,RR} = Multiplication(f,space(g))*g




## Derivative

# override map definition
Derivative(S::Hardy{s,DD},k::Integer) where {s,DD<:Circle} = ConcreteDerivative(S,k)
Derivative(S::Hardy{s,DD},k::Integer) where {s,DD<:PeriodicInterval} = ConcreteDerivative(S,k)
Derivative(S::Laurent{DD,RR},k::Integer) where {DD<:Circle,RR} =
    DerivativeWrapper(InterlaceOperator(Diagonal([map(s->Derivative(s,k),S.spaces)...]),SumSpace),k)

bandinds(D::ConcreteDerivative{Hardy{s,DD,RR}}) where {s,DD<:PeriodicInterval,RR}=(0,0)
bandinds(D::ConcreteDerivative{Hardy{s,DD,RR}}) where {s,DD<:Circle,RR}=s?(0,D.order):(-D.order,0)

rangespace(D::ConcreteDerivative{S}) where {S<:Hardy}=D.space

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


getindex(D::ConcreteDerivative{Taylor{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T}  =
    T(taylor_derivative_getindex(domain(D),D.order,k,j))
getindex(D::ConcreteDerivative{Hardy{false,DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} =
    T(hardyfalse_derivative_getindex(domain(D),D.order,k,j))




## Integral



# Integral{D<:Circle}(S::Taylor{D},m)=Integral{Taylor,typeof(m),Complex{Float64}}(S,m)
#
# function Integral{D<:PeriodicInterval}(S::Hardy{false,D},m)
#     Integral{Hardy{false,D},typeof(m),Complex{Float64}}(S,m)
# end


Integral(S::Hardy{s,DD,RR},k::Integer) where {s,DD<:Circle,RR} = ConcreteIntegral(S,k)

bandinds(D::ConcreteIntegral{Taylor{DD,RR}}) where {DD<:Circle,RR} = (-D.order,0)
rangespace(Q::ConcreteIntegral{Hardy{s,DD,RR}}) where {s,DD<:Circle,RR} = Q.space

function getindex(D::ConcreteIntegral{Taylor{DD,RR}},k::Integer,j::Integer) where {DD<:Circle,RR}
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


function Integral(S::SubSpace{Hardy{false,DD,RR},UnitCount{Int64},DD,RR},k::Integer) where {DD<:Circle,RR}
    if first(S.indexes) == k+1
        ConcreteIntegral(S,k)
    else
        @assert first(S.index) > k+1
        S2=S.space|(k+1:∞)
        IntegralWrapper(ConcreteIntegral(S2,k)*Converion(S,S2),k)
    end
end

bandinds(D::ConcreteIntegral{SubSpace{Hardy{false,DD,RR},UnitCount{Int64},DD,RR}}) where {DD<:Circle,RR} =
    (0,0)

rangespace(D::ConcreteIntegral{SubSpace{Hardy{false,DD,RR},UnitCount{Int64},DD,RR}}) where {DD<:Circle,RR} =
    D.space.space

function getindex(D::ConcreteIntegral{SubSpace{Hardy{false,DD,RR},UnitCount{Int64},DD,RR},OT,TT},
                 k::Integer,j::Integer) where {RR,OT,TT,DD<:Circle}
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



bandinds(D::ConcreteIntegral{Hardy{false,DD,RR}}) where {DD<:PeriodicInterval,RR} = (0,0)
rangespace(D::ConcreteIntegral{Taylor{DD,RR}}) where {DD<:PeriodicInterval,RR} = D.space


function getindex(D::ConcreteIntegral{Hardy{false,DD,RR}},k::Integer,j::Integer) where {DD<:PeriodicInterval,RR}
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



bandinds(D::ConcreteIntegral{SubSpace{Taylor{DD,RR},UnitCount{Int64},DD,RR}}) where {RR,DD<:PeriodicInterval} =
    (0,0)
rangespace(D::ConcreteIntegral{SubSpace{Taylor{DD,RR},UnitCount{Int64},DD,RR}}) where {RR,DD<:PeriodicInterval} =
    D.space

function getindex(D::ConcreteIntegral{SubSpace{Taylor{DD,RR},UnitCount{Int64},DD,RR}},
                 k::Integer,j::Integer) where {RR,DD<:PeriodicInterval}
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
        DefiniteIntegral(S::$SP{D,R}) where {D,R} =
            ConcreteDefiniteIntegral{typeof(S),prectype(S)}(S)
        DefiniteLineIntegral(S::$SP{D,R}) where {D,R} =
            ConcreteDefiniteLineIntegral{typeof(S),real(prectype(S))}(S)
    end
end

getindex(Σ::ConcreteDefiniteIntegral{Taylor{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Hardy{false,D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Laurent{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Taylor{D,R},T},k::Integer) where {T,D<:Circle,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Hardy{false,D,R},T},k::Integer) where {T,D<:Circle,R} =
    k == 1? T(complexlength(domain(Σ))) :zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Laurent{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k == 2? T(complexlength(domain(Σ))) :zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{Taylor{D,R},T},k::Integer) where {T,D,R} =
    k == 1? T(arclength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{Hardy{false,D,R},T},k::Integer) where {T,D,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{Laurent{D,R},T},k::Integer) where {T,D,R} =
    k == 1? T(arclength(domain(Σ))) : zero(T)

bandinds(Σ::ConcreteDefiniteIntegral{Taylor{D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Hardy{false,D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Laurent{D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Taylor{D,R}}) where {D<:Circle,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Hardy{false,D,R}}) where {D<:Circle,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Laurent{D,R}}) where {D<:Circle,R} = 0,1
bandinds(Σ::ConcreteDefiniteLineIntegral{Taylor{D,R}}) where {D,R} = 0,0
bandinds(Σ::ConcreteDefiniteLineIntegral{Hardy{false,D,R}}) where {D,R} = 0,0
bandinds(Σ::ConcreteDefiniteLineIntegral{Laurent{D,R}}) where {D,R} = 0,0


## reverse orientation

conversion_type(A::Laurent{DD,RR},B::Laurent{DD,RR}) where {DD<:Circle,RR}=domain(A).orientation?A:B
function Conversion(A::Laurent{DD,RR},B::Laurent{DD,RR}) where {DD,RR}
    @assert domain(A) == reverse(domain(B))
    ConversionWrapper(SpaceOperator(
        InterlaceOperator(Diagonal([eye(1),PermutationOperator([2,1])]))
    ,A,B))
end
