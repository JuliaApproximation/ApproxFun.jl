# DirichletSpaces

export ChebyshevDirichlet

#TODO: Type of Segment
struct ChebyshevDirichlet{left,right,D,R} <: PolynomialSpace{D,R}
    domain::D
    ChebyshevDirichlet{left,right,D,R}(d) where {left,right,D,R} = new(d)
    ChebyshevDirichlet{left,right,D,R}() where {left,right,D,R} = new(D())
end

for TYP in (:Number,:AbstractArray,:Vec,:Fun)
    @eval evaluate(f::AbstractVector,S::ChebyshevDirichlet,x::$TYP) =
        evaluate(Fun(Fun(S,f),canonicalspace(S)),x)
end
ChebyshevDirichlet{l,r}() where {l,r} =
    ChebyshevDirichlet{l,r,Segment{Float64},Float64}()
ChebyshevDirichlet{l,r}(d::Domain) where {l,r} =
    ChebyshevDirichlet{l,r,typeof(d),real(prectype(d))}(d)

spacescompatible{l,r,D,R}(a::ChebyshevDirichlet{l,r,D,R},b::ChebyshevDirichlet{l,r,D,R}) =
    domainscompatible(a,b)

ChebyshevDirichlet() = ChebyshevDirichlet{1,1,Segment{Float64},Float64}()
ZeroChebyshevDirichlet(d) =
    ChebyshevDirichlet{1,1,Segment{Float64}}(d)|(3:∞)
ZeroChebyshevDirichlet() =
    ChebyshevDirichlet{1,1,Segment{Float64}}()|(3:∞)

canonicalspace(S::ChebyshevDirichlet) = Chebyshev(domain(S))

setdomain{l,r}(S::ChebyshevDirichlet{l,r},d::Domain) = ChebyshevDirichlet{l,r}(d)


# These are used to make sure Chebyshev comes first

Base.isless(a::Chebyshev,b::ChebyshevDirichlet)=true
<(a::Chebyshev,b::ChebyshevDirichlet)=true
<=(a::Chebyshev,b::ChebyshevDirichlet)=true
>(a::Chebyshev,b::ChebyshevDirichlet)=false
>=(a::Chebyshev,b::ChebyshevDirichlet)=false

Base.isless(a::ChebyshevDirichlet,b::Chebyshev)=false
<(a::ChebyshevDirichlet,b::Chebyshev)=false
<=(a::ChebyshevDirichlet,b::Chebyshev)=false
>(a::ChebyshevDirichlet,b::Chebyshev)=true
>=(a::ChebyshevDirichlet,b::Chebyshev)=true

## coefficients

# converts to f_0 T_0 + f_1 T_1 +  \sum f_k (T_k - T_{k-2})

function dirichlettransform!(w::AbstractVector)
    for k=length(w)-2:-1:1
        @inbounds w[k] += w[k+2]
    end

    w
end

function idirichlettransform!(w::AbstractVector)
    for k=3:length(w)
        @inbounds w[k-2]-= w[k]
    end

    w
end


# converts to f_0 T_0 + \sum f_k (T_k ± T_{k-1})

function idirichlettransform!(s::Bool,w::AbstractVector)
    for k=2:length(w)
        @inbounds w[k-1]+= (s?-1:1)*w[k]
    end

    w
end


function dirichlettransform!(s::Bool,w::AbstractVector)
    for k=length(w)-1:-1:1
        @inbounds w[k] += (s?1:-1)*w[k+1]
    end

    w
end



coefficients(v::AbstractVector,::Chebyshev,::ChebyshevDirichlet{1,1})=dirichlettransform!(copy(v))
coefficients(v::AbstractVector,::Chebyshev,::ChebyshevDirichlet{0,1})=dirichlettransform!(true,copy(v))
coefficients(v::AbstractVector,::Chebyshev,::ChebyshevDirichlet{1,0})=dirichlettransform!(false,copy(v))

coefficients(v::AbstractVector,::ChebyshevDirichlet{1,1},::Chebyshev)=idirichlettransform!(copy(v))
coefficients(v::AbstractVector,::ChebyshevDirichlet{0,1},::Chebyshev)=idirichlettransform!(true,copy(v))
coefficients(v::AbstractVector,::ChebyshevDirichlet{1,0},::Chebyshev)=idirichlettransform!(false,copy(v))

# recurrence
recα{T}(::Type{T},::ChebyshevDirichlet{1,1},n) = zero(T)
recβ{T}(::Type{T},::ChebyshevDirichlet{1,1},n) = n==1 ? one(T) : one(T)/2
recγ{T}(::Type{T},::ChebyshevDirichlet{1,1},n) = n==2 ? one(T) : (n==3 ? zero(T) : one(T)/2)

## Dirichlet Conversion operators

Conversion(D::ChebyshevDirichlet,C::Chebyshev)=ConcreteConversion(D,C)

getindex{D,R,CC<:Chebyshev,T}(C::ConcreteConversion{ChebyshevDirichlet{1,0,D,R},CC,T},k::Integer,j::Integer) =
    j==k || j==k+1 ? one(T) : zero(T)

getindex{D,R,CC<:Chebyshev,T}(C::ConcreteConversion{ChebyshevDirichlet{0,1,D,R},CC,T},k::Integer,j::Integer) =
    j==k ? one(T) : ( j==k+1? -one(T) : zero(eltype(C)))

getindex{D,R,CC<:Chebyshev,T}(C::ConcreteConversion{ChebyshevDirichlet{1,1,D,R},CC,T},k::Integer,j::Integer) =
    j==k ? one(T) : ( j==k+2? -one(T) : zero(eltype(C)))

function getindex{D,R,CC<:Chebyshev,T}(C::ConcreteConversion{ChebyshevDirichlet{2,2,D,R},CC,T},k::Integer,j::Integer)
    if j==k
        one(T)
    elseif j==k+4
        2one(T)*(k+1)/k-1
    elseif k≥3 && j==k+2
        -2one(T)*(k-1)/(k-2)
    else
        zero(T)
    end
end


bandinds{D,R,C<:Chebyshev}(::ConcreteConversion{ChebyshevDirichlet{1,0,D,R},C})=0,1
bandinds{D,R,C<:Chebyshev}(::ConcreteConversion{ChebyshevDirichlet{0,1,D,R},C})=0,1
bandinds{D,R,C<:Chebyshev}(::ConcreteConversion{ChebyshevDirichlet{1,1,D,R},C})=0,2
bandinds{D,R,C<:Chebyshev}(::ConcreteConversion{ChebyshevDirichlet{2,2,D,R},C})=0,4

conversion_rule(b::ChebyshevDirichlet,a::Chebyshev)=b

# return the space that has banded Conversion to the other
# function conversion_rule(a::ChebyshevDirichlet,b::Ultraspherical)
#     @assert domainscompatible(a,b)
#
#     a
# end




## Evaluation Functional


bandinds{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{1,0,D,R},Bool}) = 0,B.x?∞:0
bandinds{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{0,1,D,R},Bool}) = 0,B.x?0:∞
bandinds{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{1,1,D,R},Bool}) = 0,1

function getindex{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{1,0,D,R},Bool},kr::Range)
    d = domain(B)

    if B.x == false && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    elseif B.x == true && B.order == 0
        Float64[k==1?1.0:2.0 for k=kr]
    else
        (Evaluation(d,B.x,B.order)*Conversion(domainspace(B)))[kr]
    end
end

function getindex{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{0,1,D,R},Bool},kr::Range)
    S = Space(domain(B))


    if B.x == true && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    elseif B.x == false && B.order == 0
        Float64[k==1?1.0:-(-1)^k*2.0 for k=kr]
    else
        (Evaluation(S,B.x,B.order)*Conversion(domainspace(B),S))[kr]
    end
end

function getindex{D,R}(B::ConcreteEvaluation{ChebyshevDirichlet{1,1,D,R},Bool},kr::Range)
    S = Space(domain(B))

    if B.x == false && B.order == 0
        Float64[k==1?1.0:(k==2?-1.0:0.0) for k=kr]
    elseif B.x == true && B.order == 0
        Float64[k==1||k==2?1.0:0.0 for k=kr]
    else
        getindex(Evaluation(S,B.x,B.order)*Conversion(domainspace(B),S),kr)
    end
end

function Evaluation(sp::ChebyshevDirichlet,x::Float64,ord::Integer)
    S=Space(domain(sp))
    EvaluationWrapper(sp,x,ord,Evaluation(S,x,ord)*Conversion(sp,S))
end
