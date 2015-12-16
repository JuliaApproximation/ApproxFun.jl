


ToeplitzOperator{DD}(f::Fun{Laurent{DD}})=ToeplitzOperator(f.coefficients[2:2:end],
                                                    f.coefficients[1:2:end])
LaurentOperator{DD}(f::Fun{Laurent{DD}})=LaurentOperator(f.coefficients[3:2:end],
                                                    f.coefficients[[1;2:2:end]])



##Taylor


bandinds{DD}(M::Multiplication{Taylor{DD},Taylor{DD}})=1-length(M.f),0
rangespace{DD}(M::Multiplication{Taylor{DD},Taylor{DD}})=domainspace(M)
addentries!{DD}(M::Multiplication{Taylor{DD},Taylor{DD}},A,k,::Colon)=addentries!(ToeplitzOperator(M.f.coefficients[2:end],[M.f.coefficients[1]]),A,k,:)


## Evaluation

getindex{DD}(T::Evaluation{Taylor{DD},Complex{Float64},Complex{Float64}},cols::Range)=mappoint(domain(T),Circle(),T.x).^(cols-1)


## Multiplication

Multiplication{DD}(f::Fun{Laurent{DD}},sp::Laurent{DD})=defaultMultiplication(f,sp)
bandinds{DD}(M::Multiplication{Laurent{DD},Laurent{DD}})=bandinds(LaurentOperator(M.f))
rangespace{DD}(M::Multiplication{Laurent{DD},Laurent{DD}})=domainspace(M)
addentries!{DD}(M::Multiplication{Laurent{DD},Laurent{DD}},A,k,::Colon)=addentries!(LaurentOperator(M.f),A,k,:)





## Derivative

# override map definition
Derivative{s,DD<:Circle}(S::Hardy{s,DD},k::Integer)=Derivative{typeof(S),typeof(k),promote_type(eltype(S),eltype(DD))}(S,k)

bandinds{s,DD<:PeriodicInterval}(D::Derivative{Hardy{s,DD}})=(0,0)
bandinds{s,DD<:Circle}(D::Derivative{Hardy{s,DD}})=s?(0,D.order):(-D.order,0)

rangespace{S<:Hardy}(D::Derivative{S})=D.space

function taylor_derivative_addentries!(d::PeriodicInterval,m,A,kr::Range)
    C=2/(d.b-d.a)*π*im
    for k=kr
        A[k,k] += (C*(k-1))^m
    end
    A
end

function hardyfalse_derivative_addentries!(d::PeriodicInterval,m,A,kr::Range)
    C=2/(d.b-d.a)*π*im
    for k=kr
        A[k,k] += (-C*k)^m
    end
    A
end



function taylor_derivative_addentries!(d::Circle,m,A,kr::Range)
    C=d.radius^(-m)

    for k=kr
        D=k
        for j=k+1:k+m-1
          D*=j
        end
        A[k,k+m] += C*D
    end

    A
end

function hardyfalse_derivative_addentries!(d::Circle,m,A,kr::Range)
    C=(-d.radius)^(-m)

    for k=max(m+1,kr[1]):kr[end]
        D=k-m
        for j=k-m+1:k-1
          D*=j
        end
        A[k,k-m] += C*D
    end

    A
end


addentries!{DD}(D::Derivative{Taylor{DD}},A,kr::Range,::Colon)=taylor_derivative_addentries!(domain(D),D.order,A,kr)
addentries!{DD}(D::Derivative{Hardy{false,DD}},A,kr::Range,::Colon)=hardyfalse_derivative_addentries!(domain(D),D.order,A,kr)




## Integral



# Integral{D<:Circle}(S::Taylor{D},m)=Integral{Taylor,typeof(m),Complex{Float64}}(S,m)
#
# function Integral{D<:PeriodicInterval}(S::Hardy{false,D},m)
#     Integral{Hardy{false,D},typeof(m),Complex{Float64}}(S,m)
# end


Integral{s,DD<:Circle}(S::Hardy{s,DD},k::Integer)=Integral{typeof(S),typeof(k),promote_type(eltype(S),eltype(DD))}(S,k)

bandinds{DD<:Circle}(D::Integral{Taylor{DD}})=(-D.order,0)
rangespace{s,DD<:Circle}(Q::Integral{Hardy{s,DD}})=Q.space

function addentries!{DD<:Circle}(D::Integral{Taylor{DD}},A,kr::Range,::Colon)
    d=domain(D)
    m=D.order

    C=d.radius^m

    for k=max(m+1,kr[1]):kr[end]
        D=k-m
        for j=k-m+1:k-1
          D*=j
        end
        A[k,k-m] += C/D
    end

    A
end


Integral{n,T,DD<:Circle}(S::SliceSpace{n,1,Hardy{false,DD},T,DD,1},k::Integer)=Integral{typeof(S),typeof(k),promote_type(eltype(S),eltype(DD))}(S,k)

function bandinds{n,T,DD<:Circle}(D::Integral{SliceSpace{n,1,Hardy{false,DD},T,DD,1}})
    @assert D.order==n
    (0,0)
end
rangespace{n,T,DD<:Circle}(D::Integral{SliceSpace{n,1,Hardy{false,DD},T,DD,1}})=D.space.space

function addentries!{n,T,DD<:Circle}(D::Integral{SliceSpace{n,1,Hardy{false,DD},T,DD,1}},A,kr::Range,::Colon)
    d=domain(D)
    m=D.order

    C=(-d.radius)^m

    for k=kr
        D=k
        for j=k+1:k+m-1
          D*=j
        end
        A[k,k] += C/D
    end

    A
end



bandinds{DD<:PeriodicInterval}(D::Integral{Hardy{false,DD}})=(0,0)
rangespace{DD<:PeriodicInterval}(D::Integral{Taylor{DD}})=D.space


function addentries!{DD<:PeriodicInterval}(D::Integral{Hardy{false,DD}},A,kr::Range,::Colon)
    d=domain(D)
    m=D.order

    C=2/(d.b-d.a)*π*im
    for k=kr
        A[k,k] += (-C*k)^(-m)
    end
    A
end



bandinds{n,T,DD<:PeriodicInterval}(D::Integral{SliceSpace{n,1,Taylor{DD},T,DD,1}})=(0,0)
rangespace{n,T,DD<:PeriodicInterval}(D::Integral{SliceSpace{n,1,Taylor{DD},T,DD,1}})=D.space

function addentries!{n,T,DD<:PeriodicInterval}(D::Integral{SliceSpace{n,1,Taylor{DD},T,DD,1}},A,kr::Range,::Colon)
    d=domain(D)
    m=D.order

    C=2/(d.b-d.a)*π*im
    for k=kr
        A[k,k] += (C*(k+n-1))^(-m)
    end
    A
end



## Definite integral

for TYP in (:DefiniteIntegral,:DefiniteLineIntegral)
    @eval $TYP{DD<:Circle}(S::Laurent{DD})=$TYP{typeof(S),promote_type(eltype(S),eltype(DD))}(S)
end

function getindex{T,DD<:PeriodicInterval}(Σ::DefiniteIntegral{Laurent{DD},T},kr::Range)
    d = domain(Σ)
    T[k == 1?  d.b-d.a : zero(T) for k=kr]
end

getindex{T,DD<:Circle}(Σ::DefiniteIntegral{Laurent{DD},T},kr::Range)=T[k == 2?  2domain(Σ).radius*π*im :zero(T) for k=kr]

datalength{DD<:PeriodicInterval}(Σ::DefiniteIntegral{Laurent{DD}})=1
datalength{DD<:Circle}(Σ::DefiniteIntegral{Laurent{DD}})=2


function getindex{T,DD<:PeriodicInterval}(Σ::DefiniteLineIntegral{Laurent{DD},T},kr::Range)
    d = domain(Σ)
    T[k == 1?  d.b-d.a : zero(T) for k=kr]
end

getindex{T,DD<:Circle}(Σ::DefiniteLineIntegral{Laurent{DD},T},kr::Range)=T[k == 1?  2domain(Σ).radius*π : zero(T) for k=kr]
datalength{DD<:PeriodicInterval}(Σ::DefiniteLineIntegral{Laurent{DD}})=1
datalength{DD<:Circle}(Σ::DefiniteLineIntegral{Laurent{DD}})=2




## reverse orientation

conversion_type{DD<:Circle}(A::Laurent{DD},B::Laurent{DD})=domain(A).orientation?A:B
function Conversion{DD}(A::Laurent{DD},B::Laurent{DD})
    @assert domain(A) == reverse(domain(B))
    ConversionWrapper(SpaceOperator(
        BlockOperator(eye(1),PermutationOperator([2,1]))
    ,A,B))
end
