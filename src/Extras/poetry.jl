#####
# This includes short-hands that are convenient but perhaps confusing
#####


export chebyshevt,chebyshevu,legendre,∫,⨜,⨍,∇,Δ,ChebyshevWeight,𝕀,ℝ,𝕌,𝒟

## Constructors

Fun()=Fun(identity)
Fun(d::Domain)=Fun(identity,d)
Fun(d::Space)=Fun(identity,d)

## Chebyshev & Legendre polynomials

chebyshevt{T<:Number}(n::Int,d::Interval{T}) = Fun(Chebyshev(d),[zeros(T,n);one(T)])
chebyshevu{T<:Number}(n::Int,d::Interval{T}) =
    mod(n,2) == 1 ? Fun(Chebyshev(d),interlace(zeros(T,div(n+2,2)),2ones(T,div(n+2,2)))) :
                    Fun(Chebyshev(d),interlace(2ones(T,div(n+2,2)),zeros(T,div(n+2,2)))[1:n+1]-[one(T);zeros(T,n)])
legendre{T<:Number}(n::Int,d::Interval{T}) = Fun(Legendre(d),[zeros(T,n);one(T)])

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



# use conj(f.') for ArraySpace
Base.ctranspose(f::Fun)=differentiate(f)


∫(f::Fun)=integrate(f)
⨜(f::Fun)=cumsum(f)

for OP in (:Σ,:∮,:⨍,:⨎)
    @eval $OP(f::Fun)=sum(f)
end

∇(F::MultivariateFun) = grad(F)
Base.dot{M<:MultivariateFun}(∇::Function,F::Vector{M}) = div(F)
Base.cross{M<:MultivariateFun}(∇::Function,F::Vector{M}) = curl(F)


## Domains

const 𝕀 = Interval()
const ℝ = Line()
const 𝕌 = Circle()

𝒟 = Derivative()
Δ = Laplacian()
