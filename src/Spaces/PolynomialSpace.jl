
## Orthogonal polynomials

abstract PolynomialSpace{D} <: RealUnivariateSpace{D}



Multiplication{U<:PolynomialSpace}(f::Fun{U},sp::PolynomialSpace)=ConcreteMultiplication(f,sp)
bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V})=(1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V})=domainspace(M)


# All polynomials contain constant
union_rule(A::ConstantSpace,B::PolynomialSpace)=B
Base.promote_rule{T<:Number,S<:PolynomialSpace,V}(::Type{Fun{S,V}},::Type{T})=Fun{S,promote_type(V,T)}
Base.promote_rule{T<:Number,S<:PolynomialSpace}(::Type{Fun{S}},::Type{T})=Fun{S,T}

## Evaluation

function evaluate(f::AbstractVector,S::PolynomialSpace,x)
    if x in domain(S)
        clenshaw(S,f,tocanonical(S,x))
    else
        zero(eltype(f))
    end
end

evaluate(f::AbstractVector,S::PolynomialSpace,x::AbstractArray)=map(y->evaluate(f,S,y),x)

# we need the ... for multi-dimensional
evaluate(f::AbstractVector,S::PolynomialSpace,x,y,z...)=
    evaluate(f,S,Vec(x,y,z...))

#evaluate(f::AbstractVector,S::PolynomialSpace,x::AbstractArray)=map(y->evaluate(f,S,y),x)

function evaluate(f::AbstractVector,S::PolynomialSpace,x::Fun)
    if issubset(Interval(minimum(x),maximum(x)),domain(S))
        clenshaw(S,f,tocanonical(S,x))
    else
        error("Implement splitatpoints for evaluate ")
    end
end



######
# Recurrence encodes the recurrence coefficients
# or equivalently multiplication by x
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

function getindex{S,T}(R::Recurrence{S,T},k::Integer,j::Integer)
    if j==k-1
        recβ(T,R.space,k-1)
    elseif j==k
        recα(T,R.space,k)
    elseif j==k+1
        recγ(T,R.space,k+1)
    else
        zero(T)
    end
end

######
# JacobiZ encodes [BasisFunctional(1);(J-z*I)[2:end,:]]
# where J is the Jacobi operator
######
immutable JacobiZ{S,T} <: TridiagonalOperator{T}
    space::S
    z::T
end

JacobiZ(sp,z)=(T = promote_type(eltype(sp),eltype(domain(sp)),typeof(z)); JacobiZ{typeof(sp),T}(sp,T(z)))

Base.convert{T,S}(::Type{BandedOperator{T}},J::JacobiZ{S})=JacobiZ{S,T}(J.space,J.z)


#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####

function getindex{S,T}(J::JacobiZ{S,T},k::Integer,j::Integer)
    if j==k-1
        recγ(T,J.space,k)
    elseif j==k
        k == 1 ? one(T) : recα(T,J.space,k)-J.z
    elseif j==k+1 && k > 1
        recβ(T,J.space,k)
    else
        zero(T)
    end
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


getindex{PS<:PolynomialSpace,T,C<:PolynomialSpace}(M::ConcreteMultiplication{C,PS,T},k::Integer,j::Integer) = M[k:k,j:j][1,1]


function Base.copy{PS<:PolynomialSpace,V,T,C<:PolynomialSpace}(S::SubBandedMatrix{T,ConcreteMultiplication{C,PS,V,T},
                                                                            Tuple{UnitRange{Int},UnitRange{Int}}})
    M=parent(S)
    kr,jr=parentindexes(S)

    A=bzeros(S)

    a=coefficients(M.f)

    shft=bandshift(A)

    for k=kr ∩ jr
        A[k-kr[1]+1,k-jr[1]+1]=a[1]
    end

    if length(a) > 1
        sp=M.space
        fsp=space(M.f)
        jkr=max(1,min(kr[1],jr[1])-length(a)+1):max(kr[end],jr[end])+length(a)-1

        #Multiplication is transpose
        J=Recurrence(sp)[jkr,jkr]


        # the sub ranges of jkr that correspond to kr, jr
        kr2,jr2=kr-jkr[1]+1,jr-jkr[1]+1


        C0=beye(size(J,1),size(J,2),0,0)
        C1=(1/recβ(T,fsp,1))*J-(recα(T,fsp,1)/recβ(T,fsp,1))*C0

        BLAS.axpy!(a[2],@compat view(C1,kr2,jr2),A)

        for k=1:length(a)-2
            rβ=recβ(T,fsp,k+1)
            C1,C0=(1/rβ)*J*C1-(recα(T,fsp,k+1)/rβ)*C1-(recγ(T,fsp,k+1)/rβ)*C0,C1
            BLAS.axpy!(a[k+2],@compat view(C1,kr2,jr2),A)
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
