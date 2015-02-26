
## Orthogonal polynomials

abstract PolynomialSpace <: IntervalSpace

bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=domainspace(M)


######
# Recurrence encodes the recurrence coefficients
# or equivalentally multiplication by x
######
immutable Recurrence{S,T} <: TridiagonalOperator{T}
    space::S
end 

Recurrence(sp)=Recurrence{typeof(sp),Float64}(sp)

Base.convert{T,S}(::Type{BandedOperator{T}},J::Recurrence{S})=Recurrence{S,T}(J.space)


function addentries!{S,T}(R::Recurrence{S,T},A,kr::Range)
    for k=kr
        A[k,k-1]=recβ(T,R.space,k-1)
        A[k,k]  =recα(T,R.space,k)        
        A[k,k+1]=recγ(T,R.space,k+1)        
    end
    A
end