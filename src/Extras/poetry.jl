#####
# This includes short-hands that are convenient but perhaps confusing
#####


export chebyshevt,chebyshevu,legendre,∫,⨜,⨍,∇,Δ,ChebyshevWeight

## Constructors

Fun()=Fun(identity)
Fun(d::Domain)=Fun(identity,d)
Fun(d::Space)=Fun(identity,d)

## Chebyshev & Legendre polynomials

chebyshevt{T<:Number}(n::Int,d::Interval{T}) = Fun([zeros(T,n),one(T)],Chebyshev(d))
chebyshevu{T<:Number}(n::Int,d::Interval{T}) = mod(n,2) == 1 ? Fun(interlace(zeros(T,div(n+2,2)),2ones(T,div(n+2,2))),Chebyshev(d)) : Fun(interlace(2ones(T,div(n+2,2)),zeros(T,div(n+2,2)))[1:n+1]-[one(T),zeros(T,n)],Chebyshev(d))
legendre{T<:Number}(n::Int,d::Interval{T}) = Fun([zeros(T,n),one(T)],Legendre(d))

for poly in (:chebyshevt,:chebyshevu,:legendre)
    @eval begin
        $poly{T<:Number}(n::Int,a::T,b::T) = $poly(n,Interval(a,b))
        $poly{T<:Number}(::Type{T},n::Int) = $poly(n,Interval{T}())
        $poly(n::Int) = $poly(Float64,n)
        $poly{T<:Number}(n::Range,d::Interval{T}) = map(i->$poly(i,d),n)
    end
end

ChebyshevWeight(d,k)=k==0?JacobiWeight(-0.5,-0.5,d):JacobiWeight(0.5,0.5,d)
ChebyshevWeight(d)=ChebyshevWeight(d,0)
ChebyshevWeight(k::Integer)=ChebyshevWeight(Interval(),k)
ChebyshevWeight()=ChebyshevWeight(0)

# shorthand for second order

ivp(d) = ivp(d,2)
bvp(d) = bvp(d,2)

## diff

# diff(::Fun{ArraySpace}) is left as array diff

Base.diff{S,T}(f::Fun{S,T},n...)=differentiate(f,n...)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)

Base.diff(d::Space,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

Base.diff(d::Union{ProductDomain,TensorSpace},k::Integer)=Derivative(d,k)

# use conj(f.') for ArraySpace
Base.ctranspose(f::Fun)=differentiate(f)


∫(f::Fun)=integrate(f)
⨜(f::Fun)=cumsum(f)

for OP in (:Σ,:∮,:⨍,:⨎)
    @eval $OP(f::Fun)=sum(f)
end

∇(F::MultivariateFun) = grad(F)
Δ(F::MultivariateFun) = lap(F)
Base.dot{M<:MultivariateFun}(∇::Function,F::Vector{M}) = div(F)
Base.cross{M<:MultivariateFun}(∇::Function,F::Vector{M}) = curl(F)
