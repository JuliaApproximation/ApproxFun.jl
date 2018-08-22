## CurveSpace


Space(d::PeriodicCurve{S}) where {S<:Fourier}=Fourier(d)
Space(d::PeriodicCurve{S}) where {S<:Laurent}=Laurent(d)

#TODO: Make type stable
Curve(f::Fun) = isa(domain(f),IntervalDomain) ? IntervalCurve(f) : PeriodicCurve(f)

# function evaluate{C<:Curve,TT}(f::AbstractVector,S::Space{TT,C},x::Number)
#     rts=roots(domain(S).curve-x)
#     @assert length(rts)==1
#     evaluate(Fun(f,setdomain(S,canonicaldomain(S))),first(rts))
# end



Fun(::typeof(identity), d::Space{<:Curve}) = Fun(setdomain(space(domain(d).curve),domain(d)),
                                                domain(d).curve.coefficients)

# Bernstein polynomials are given by:
#
# b_{ν,n}(x) = binomial(n,ν)/2^n*(1+x)^ν(1-x)^{n-v},  for x ∈ [-1,1], ν = 0,…,n.
#
export Bernstein, Bézier

struct Bernstein{order,T,R} <: Space{T,R}
    domain::Segment{T}
    Bernstein{order,T,R}(d) where {order,T,R} = new(d)
    Bernstein{order,T,R}() where {order,T,R} = new(Segment{T}())
end

const Bézier = Bernstein # option+e e gives é

Bernstein{O}() where {O} = Bernstein{O,Float64,Float64}()
Bernstein{O}(d::Domain) where {O} = Bernstein{O,eltype(d),real(prectype(d))}(d)

order(::Bernstein{O}) where {O} = O
order(::Type{Bernstein{O,T,R}}) where {O,T,R} = O
dimension(::Bernstein{O}) where {O} = O+1
dimension(::Type{Bernstein{O,T,R}}) where {O,T,R} = O+1

canonicalspace(B::Bernstein) = Chebyshev(domain(B))
canonicaldomain(B::Bernstein{O,T}) where {O,T} = Segment{T}()

spacescompatible(a::Bernstein{O,T},b::Bernstein{O,T}) where {O,T}=domainscompatible(a,b)

setdomain(S::Bernstein{O},d::Domain) where {O}=Bernstein{O}(d)

Fun(::typeof(identity), B::Bernstein{order,T}) where {order,T} =
    Fun(B,collect(-one(T):2one(T)/order:one(T)))

evaluate(f::AbstractVector,S::Bernstein,z) = decasteljau(f,S,tocanonical(S,z))

# de Casteljau's numerically stable evaluation of Bernstein polynomials

function decasteljau(f::AbstractVector,S::Bernstein,z)
    β,omz,opz = copy(f),one(z)-z,one(z)+z
    for i=order(S):-1:1
        for j=1:i
            β[j] = (β[j]*omz+β[j+1]*opz)/2
        end
    end
    β[1]
end

function decasteljau(f::AbstractVector,S::Bernstein,z::AbstractVector)
    βinit = copy(f)
    β = copy(βinit)
    T = promote_type(eltype(f),eltype(z))
    ret = zeros(T,length(z))
    for k=1:length(z)
        omzk,opzk = one(z[k])-z[k],one(z[k])+z[k]
        for i=order(S):-1:1
            for j=1:i
                β[j] = (β[j]*omzk+β[j+1]*opzk)/2
            end
        end
        ret[k] = β[1]
        copyto!(β,βinit)
    end
    ret
end

splitbernstein(f::Fun{B,T},z) where {B<:Bernstein,T} = splitbernstein(coefficients(f),space(f),z)

function splitbernstein(f::AbstractVector,S::Bernstein,z)
    β,omz,opz = copy(f),one(z)-z,one(z)+z
    β1,β2 = copy(β),copy(β)
    for i=order(S):-1:1
        for j=1:i
            β[j] = (β[j]*omz+β[j+1]*opz)/2
        end
        β1[order(S)+2-i] = β[1]
        β2[i] = β[i]
    end
    Fun(PiecewiseSpace(Bernstein{order(S)}(Segment(first(domain(S)),z)),Bernstein{order(S)}(Segment(z,last(domain(S))))),
        interlace(β1,β2))
end

Fun(f::Function,S::Bernstein) = Fun(Fun(f,canonicalspace(S),dimension(S)),S)

# Connection coefficients with Chebyshev polynomials are derived in
#
#    A. Rababah, "Transformation of Chebyshev–Bernstein polynomial basis", Comp. Meth. Appl. Math. 3:608–622, 2003.
#
function coefficients(f::AbstractVector,a::Chebyshev,b::Bernstein)
    if domain(a) == domain(b)
        n = length(f)-1
        @assert n == order(b)
        M = zeros(eltype(f),n+1,n+1)
        for j=0:n
            for k=0:n
                temp = zero(eltype(f))
                for i=max(0,j+k-n):min(j,k)
                    temp += (-1)^(k-i)*binomial(2k,2i)*binomial(n-k,j-i)
                end
                M[j+1,k+1] = temp/binomial(n,j)
            end
        end
        M*f
    else
        defaultcoefficients(f,a,b)
    end
end

function coefficients(f::AbstractVector,a::Bernstein,b::Chebyshev)
    if domain(a) == domain(b)
        n = length(f)-1
        @assert n == order(a)
        Minv = zeros(eltype(f),n+1,n+1)
        for j=0:n
            for k=0:n
                temp = zero(eltype(f))
                for i=0:j
                    temp += (-1)^(j-i)*binomial(2j,2i)*binomial(2k+2i,k+i)*binomial(2n+2j-2k-2i,n+j-k-i)/binomial(n+j,k+i)
                end
                Minv[j+1,k+1] = (2-FastTransforms.δ(j,0))/4^(n+j)*binomial(n,k)*temp
            end
        end
        Minv*f
    else
        defaultcoefficients(f,a,b)
    end
end
