## This makes implementing operators simpler
# but overrided copy directly is likely to be faster


abstract TridiagonalOperator{T} <: Operator{T}
abstract SymTridiagonalOperator{T} <: TridiagonalOperator{T}
abstract BidiagonalOperator{T} <: TridiagonalOperator{T}
abstract DiagonalOperator{T} <: BidiagonalOperator{T}
# override getindex

bandinds(::TridiagonalOperator)=-1,1
bandinds(::BidiagonalOperator)=0,1
bandinds(::DiagonalOperator)=0,0



Base.Tridiagonal{T}(J::TridiagonalOperator{T},n::Integer)=Tridiagonal(T[J.γ(k) for k=2:n],
T[J.α(k) for k=1:n],
T[J.β(k) for k=1:n-1])


function symmetrize{T}(J::TridiagonalOperator{T},n::Integer)
    d=Array(T,n)
    d[1]=1
    for k=2:n
        d[k]=sqrt(J[k,k-1]/J[k-1,k])*d[k-1]
    end

   SymTridiagonal(
    T[J[k,k] for k=1:n],
    T[J[k,k+1]*d[k+1]/d[k] for k=1:n-1])
end

immutable DiagonalInverseOperator{S,T} <: DiagonalOperator{T}
    op::S
end

function DiagonalInverseOperator(B::Operator)
    @assert bandinds(B)==(0,0)
    DiagonalInverseOperator{typeof(B),promote_type(Float64,eltype(B))}(B)
end

domainspace(D::DiagonalInverseOperator)=rangespace(D.op)
rangespace(D::DiagonalInverseOperator)=domainspace(D.op)

Base.getindex(D::DiagonalInverseOperator,k::Integer,j::Integer)=k==j?1./D.op[k,k]:zero(eltype(D))


@eval Base.convert{T}(::Type{Operator{T}},F::DiagonalInverseOperator)=DiagonalInverseOperator{typeof(F.op),
                                                                                        T}(F.op)

Base.inv(B::Operator)=DiagonalInverseOperator(B)
