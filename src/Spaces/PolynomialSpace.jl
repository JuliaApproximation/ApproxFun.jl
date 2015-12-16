
## Orthogonal polynomials

abstract PolynomialSpace{D} <: RealUnivariateSpace{D}

bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V})=domainspace(M)


# All polynomials contain constant
union_rule(A::ConstantSpace,B::PolynomialSpace)=B
Base.promote_rule{T<:Number,S<:PolynomialSpace,V}(::Type{Fun{S,V}},::Type{T})=Fun{S,promote_type(V,T)}
Base.promote_rule{T<:Number,S<:PolynomialSpace}(::Type{Fun{S}},::Type{T})=Fun{S,T}

## Evaluation

evaluate(f::AbstractVector,S::PolynomialSpace,x)=clenshaw(S,f,tocanonical(S,x))

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

function addentries!{S,T}(R::Recurrence{S,T},A,kr::Range,::Colon)
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


function addentries!{US<:PolynomialSpace,PS<:PolynomialSpace,T}(M::ConcreteMultiplication{US,PS,T},A,kr::UnitRange,::Colon)
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
        addentries!(C1,a[2],A,kr,:)


        for k=1:length(a)-2
            rβ=recβ(T,fsp,k+1)
            C1,C0=(1/rβ)*J*C1-(recα(T,fsp,k+1)/rβ)*C1-(recγ(T,fsp,k+1)/rβ)*C0,C1
            addentries!(C1,a[k+2],A,kr,:)
        end
    end

    A
end




## All polynomial spaces can be converted provided spaces match

isconvertible(a::PolynomialSpace,b::PolynomialSpace)=domain(a)==domain(b)
union_rule{D}(a::PolynomialSpace{D},b::PolynomialSpace{D})=domainscompatible(a,b)?(a<b?a:b):NoSpace()   # the union of two polys is always a poly


## General clenshaw
clenshaw(sp::PolynomialSpace,c::AbstractVector,x::AbstractArray) = clenshaw(c,x,
                                            ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))
clenshaw(sp::PolynomialSpace,c::AbstractMatrix,x::AbstractArray) = clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,size(c,1),length(x)))
clenshaw(sp::PolynomialSpace,c::AbstractMatrix,x) = clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),typeof(x)),sp,size(c,1),size(c,2)))

clenshaw!(sp::PolynomialSpace,c::AbstractVector,x::AbstractVector)=clenshaw!(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(x)))

function clenshaw(sp::PolynomialSpace,c::AbstractVector,x)
    N,T = length(c),promote_type(eltype(c),typeof(x))
    TT = eltype(T)
    if isempty(c)
        return zero(T)
    end

    bk1,bk2 = zero(T),zero(T)
    A,B,C = recA(TT,sp,N-1),recB(TT,sp,N-1),recC(TT,sp,N)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[k])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[k])) # (A*x+B)*bk1+c[k]-C*bk2
        A,B,C = recA(TT,sp,k-2),recB(TT,sp,k-2),recC(TT,sp,k-1)
    end
    muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[1])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[1])) # (A*x+B)*bk1+c[1]-C*bk2
end
