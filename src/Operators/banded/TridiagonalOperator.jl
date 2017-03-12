## This makes implementing operators simpler
# but overrided BandedMatrix directly is likely to be faster


@compat abstract type TridiagonalOperator{T} <: Operator{T} end
@compat abstract type SymTridiagonalOperator{T} <: TridiagonalOperator{T} end
@compat abstract type BidiagonalOperator{T} <: TridiagonalOperator{T} end
@compat abstract type DiagonalOperator{T} <: BidiagonalOperator{T} end
# override getindex

bandinds(::TridiagonalOperator)=-1,1
bandinds(::BidiagonalOperator)=0,1
bandinds(::DiagonalOperator)=0,0



Base.Tridiagonal{T}(J::TridiagonalOperator{T},n::Integer)=Tridiagonal(T[J.γ(k) for k=2:n],
T[J.α(k) for k=1:n],
T[J.β(k) for k=1:n-1])


function symmetrize{T}(J::TridiagonalOperator{T},n::Integer)
    d=Vector{T}(n)
    d[1]=1
    for k=2:n
        d[k]=sqrt(J[k,k-1]/J[k-1,k])*d[k-1]
    end

   SymTridiagonal(
    T[J[k,k] for k=1:n],
    T[J[k,k+1]*d[k+1]/d[k] for k=1:n-1])
end


immutable DiagIteratorOperator{IT,T} <: DiagonalOperator{T}
    iterator::IT
end

DiagIteratorOperator(it) = DiagIteratorOperator{typeof(it),eltype(it)}(it)

getindex(D::DiagIteratorOperator,k::Integer,j::Integer) =
    k==j? D.iterator[k] : zero(eltype(D))

domainspace(D::DiagIteratorOperator) = ℓ⁰
rangespace(D::DiagIteratorOperator) = ℓ⁰

Base.diagm(c::AbstractCount) = DiagIteratorOperator(c)
