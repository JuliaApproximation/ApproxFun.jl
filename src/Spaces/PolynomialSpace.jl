
## Orthogonal polynomials

abstract PolynomialSpace <: IntervalSpace

bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::Multiplication{U,V})=domainspace(M)


# All polynomials contain constant
union_rule(A::ConstantSpace,B::PolynomialSpace)=B

######
# Recurrence encodes the recurrence coefficients
# or equivalentally multiplication by x
######
immutable Recurrence{S,T} <: TridiagonalOperator{T}
    space::S
end

Recurrence(sp)=Recurrence{typeof(sp),promote_type(eltype(sp),eltype(domain(sp)))}(sp)

Base.convert{T,S}(::Type{BandedOperator{T}},J::Recurrence{S})=Recurrence{S,T}(J.space)


#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####

function addentries!{S,T}(R::Recurrence{S,T},A,kr::Range)
    for k=kr
        A[k,k-1]=recβ(T,R.space,k-1)
        A[k,k]  =recα(T,R.space,k)
        A[k,k+1]=recγ(T,R.space,k+1)
    end
    A
end




#####
# Multiplication can be built from recurrence coefficients
#  multiplication by J=M[x] is Recurrence operator
#  and assume p_0 = 1, then we have
#
#   M[p_0] = 1
#   M[p_1] = (J/β_1 - α_1/β_1)*M[p_0]
#   M[p_k] = (J/β_k - α_k/β_k)*M[p_{k-1}] - γ_k/β_k*M[p_{k-2}]
#####


function addentries!{US<:PolynomialSpace,PS<:PolynomialSpace,T}(M::Multiplication{US,PS,T},A,kr::UnitRange)
    a=coefficients(M.f)
    fsp=space(M.f)
    for k=kr
        A[k,k]=a[1]
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1

        J=subview(Recurrence(domainspace(M)),jkr,jkr)
        C0=isbaeye(jkr)
        C1=(1/recβ(T,fsp,1))*J-(recα(T,fsp,1)/recβ(T,fsp,1))*C0
        addentries!(C1,a[2],A,kr)


        for k=1:length(a)-2
            rβ=recβ(T,fsp,k+1)
            C1,C0=(1/rβ)*J*C1-(recα(T,fsp,k+1)/rβ)*C1-(recγ(T,fsp,k+1)/rβ)*C0,C1
            addentries!(C1,a[k+2],A,kr)
        end
    end

    A
end




## All polynomial spaces can be converted provided spaces match

isconvertible(a::PolynomialSpace,b::PolynomialSpace)=domain(a)==domain(b)
