#####
# This includes short-hands that are convenient but perhaps confusing
#####


export ∫,⨜,⨍,chebyshevt,chebyshevu

## Constructors

Fun()=Fun(identity)
Fun(d::Domain)=Fun(identity,d)
Fun(d::FunctionSpace)=Fun(identity,d)

## Chebyshev polynomials

chebyshevt{T<:Number}(::Type{T},n::Int,a::T,b::T) = Fun([zeros(T,n),one(T)],Chebyshev([a,b]))
chebyshevu{T<:Number}(::Type{T},n::Int,a::T,b::T) = mod(n,2) == 1 ? Fun(interlace(zeros(T,div(n+2,2)),2ones(T,div(n+2,2))),Chebyshev([a,b])) : Fun(interlace(2ones(T,div(n+2,2)),zeros(T,div(n+2,2)))[1:n+1]-[one(T),zeros(T,n)],Chebyshev([a,b]))

for poly in (:chebyshevt,:chebyshevu)
    @eval begin
        $poly{T<:Number}(n::Int,a::T,b::T) = $poly(Float64,n,@compat(Float64(a)),@compat(Float64(b)))
        $poly{T<:Number}(::Type{T},n::Int) = $poly(T,n,-one(T),one(T))
        $poly(n::Int) = $poly(Float64,n)
    end
end

## diff

# diff(::Fun{ArraySpace}) is left as array diff

Base.diff{S,T}(f::Fun{S,T},n...)=differentiate(f,n...)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)

Base.diff(d::FunctionSpace,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

Base.diff(d::Union(ProductDomain,TensorSpace),k::Integer)=Derivative(d,k)

# use conj(f.') for ArraySpace
Base.ctranspose(f::Fun)=differentiate(f)


∫(f::Fun)=integrate(f)
⨜(f::Fun)=cumsum(f)

for OP in (:Σ,:∮,:⨍,:⨎)
    @eval $OP(f::Fun)=sum(f)
end
