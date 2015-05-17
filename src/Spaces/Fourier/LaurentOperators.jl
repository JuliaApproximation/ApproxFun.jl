


ToeplitzOperator(f::Fun{Laurent})=ToeplitzOperator(f.coefficients[2:2:end],
                                                    f.coefficients[1:2:end])
LaurentOperator(f::Fun{Laurent})=LaurentOperator(f.coefficients[3:2:end],
                                                    f.coefficients[[1;2:2:end]])



##Taylor


bandinds(M::Multiplication{Taylor,Taylor})=1-length(M.f),0
rangespace(M::Multiplication{Taylor,Taylor})=domainspace(M)
addentries!(M::Multiplication{Taylor,Taylor},A,k)=addentries!(ToeplitzOperator(M.f.coefficients[2:end],[M.f.coefficients[1]]),A,k)


## Evaluation

getindex(T::Evaluation{Taylor,Complex{Float64},Complex{Float64}},cols::Range)=mappoint(domain(T),Circle(),T.x).^(cols-1)


## Multiplication

bandinds(M::Multiplication{Laurent,Laurent})=bandinds(LaurentOperator(M.f))
rangespace(M::Multiplication{Laurent,Laurent})=domainspace(M)
addentries!(M::Multiplication{Laurent,Laurent},A,k)=addentries!(LaurentOperator(M.f),A,k)





## Derivative

function bandinds{s}(D::Derivative{Hardy{s}})
    d=domain(D)
    if isa(d,PeriodicInterval)
        (0,0)
    elseif isa(d,Circle)
        s?(0,D.order):(-D.order,0)
    else
        error("Derivative not defined for "*string(typeof(d)))
    end
end
rangespace{S<:Hardy}(D::Derivative{S})=D.space

function taylor_derivative_addentries!(d::PeriodicInterval,m::Integer,A,kr::Range)
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,k] += (C*(k-1))^m
    end
    A
end

function hardyfalse_derivative_addentries!(d::PeriodicInterval,m::Integer,A,kr::Range)
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,k] += (-C*k)^m
    end
    A
end



function taylor_derivative_addentries!(d::Circle,m::Integer,A,kr::Range)
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

function hardyfalse_derivative_addentries!(d::Circle,m::Integer,A,kr::Range)
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


addentries!(D::Derivative{Taylor},A,kr::Range)=taylor_derivative_addentries!(domain(D),D.order,A,kr)
addentries!(D::Derivative{Hardy{false}},A,kr::Range)=hardyfalse_derivative_addentries!(domain(D),D.order,A,kr)




## Integral



function Integral(S::Taylor,m)
    @assert isa(domain(S),Circle)
    Integral{Taylor,typeof(m),Complex{Float64}}(S,m)
end

function Integral(S::Hardy{false},m)
    @assert isa(domain(S),PeriodicInterval)
    Integral{Hardy{false},typeof(m),Complex{Float64}}(S,m)
end

function bandinds(D::Integral{Taylor})
    d=domain(D)
    @assert isa(d,Circle)
    (-D.order,0)
end
rangespace(D::Integral{Taylor})=D.space
rangespace(Q::Integral{Hardy{false}})=Q.space

function addentries!(D::Integral{Taylor},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,Circle)

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


function bandinds{n,T}(D::Integral{SliceSpace{n,1,Hardy{false},T,1}})
    d=domain(D)
    @assert isa(d,Circle)
    @assert D.order==n
    (0,0)
end
rangespace{n,T}(D::Integral{SliceSpace{n,1,Hardy{false},T,1}})=D.space.space

function addentries!{n,T}(D::Integral{SliceSpace{n,1,Hardy{false},T,1}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,Circle)

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



function bandinds(D::Integral{Hardy{false}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace(D::Integral{Taylor})=D.space


function addentries!(D::Integral{Hardy{false}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,PeriodicInterval)

    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,k] += (-C*k)^(-m)
    end
    A
end



function bandinds{n,T}(D::Integral{SliceSpace{n,1,Taylor,T,1}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace{n,T}(D::Integral{SliceSpace{n,1,Taylor,T,1}})=D.space

function addentries!{n,T}(D::Integral{SliceSpace{n,1,Taylor,T,1}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,PeriodicInterval)

    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,k] += (C*(k+n-1))^(-m)
    end
    A
end



## Definite integral

function getindex{T}(Σ::DefiniteIntegral{Laurent,T},kr::Range)
    d = domain(Σ)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 2?  2π*im*d.radius :zero(T) for k=kr]
    end
end

datalength(Σ::DefiniteIntegral{Laurent})=isa(domain(Σ),PeriodicInterval)?1:2

function getindex{T}(Σ::DefiniteLineIntegral{Laurent,T},kr::Range)
    d = domain(Σ)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 1?  2d.radius*π : zero(T) for k=kr]
    end
end

datalength(Σ::DefiniteLineIntegral{Laurent})=1

