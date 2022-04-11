#####
# This includes short-hands that are convenient but perhaps confusing
#####


export chebyshevt, chebyshevu, legendre, ∫, ⨜, ∑, ⨍, ChebyshevWeight, 𝕀, 𝕌, 𝒟, ∆

## Chebyshev & Legendre polynomials

chebyshevt(n::Int,d::IntervalOrSegment{T}) where {T<:Number} = Fun(Chebyshev(d),[zeros(T,n);one(T)])
chebyshevu(n::Int,d::IntervalOrSegment{T}) where {T<:Number} =
    mod(n,2) == 1 ? Fun(Chebyshev(d),interlace(zeros(T,div(n+2,2)),2fill(one(T),div(n+2,2)))) :
                    Fun(Chebyshev(d),interlace(2fill(one(T),div(n+2,2)),zeros(T,div(n+2,2)))[1:n+1]-[one(T);zeros(T,n)])
legendre(n::Int,d::IntervalOrSegment{T}) where {T<:Number} = Fun(Legendre(d),[zeros(T,n);one(T)])

for poly in (:chebyshevt,:chebyshevu,:legendre)
    @eval begin
        $poly(::Type{T}, n::Int) where {T<:Number} = $poly(n,ChebyshevInterval{T}())
        $poly(n::Int) = $poly(Float64,n)
        $poly(n::AbstractRange, d::IntervalOrSegment{T}) where {T<:Number} = map(i->$poly(i,d),n)
    end
end

ChebyshevWeight(d,k)=k==0 ? JacobiWeight(-0.5,-0.5,d) : JacobiWeight(0.5,0.5,d)
ChebyshevWeight(d)=ChebyshevWeight(d,0)
ChebyshevWeight(k::Integer)=ChebyshevWeight(Segment(),k)
ChebyshevWeight()=ChebyshevWeight(0)



## diff
∫(f::Fun)=integrate(f)
⨜(f::Fun)=cumsum(f)

for OP in (:∑,:∮,:⨍,:⨎)  # ∑ entered by \sum<tab>, not \Sigma<tab>
    @eval $OP(f::Fun)=sum(f)
end

∇(F::MultivariateFun) = grad(F)
dot(∇::Function,F::Vector{M}) where {M<:MultivariateFun} = div(F)
cross(∇::Function,F::Vector{M}) where {M<:MultivariateFun} = curl(F)


## Domains

const 𝕀 = ChebyshevInterval()  # \bbI<tab>
const ℝ = Line()  # \bbR<tab>
const 𝕌 = Circle()  # \bbU<tab>

𝒟 = Derivative()  # \scrD<tab>
∆ = Laplacian()  # \increment<tab>, not \Delta<tab>
